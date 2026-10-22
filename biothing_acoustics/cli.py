"""Command line: breed a desk screen for a brief and write the full design package."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import optimize
from .brief import load_brief
from .design import DeskScreen, evaluate, face_relief
from .export import write_bom_csv, write_obj, write_svg
from .materials import CORES, FELTS
from .modular import JOINTS, SUPPORTS, bill_of_materials, plan_layout
from .report import design_json, make_plots, write_report


def build_package(best, brief, out, history, shortlist):
    out.mkdir(parents=True, exist_ok=True)
    m = best.metrics
    runs = plan_layout(brief.layout, brief.module_widths)
    bom = bill_of_materials(best.screen, runs, m["surface_mass_kg_m2"], m["thickness_mm"],
                            m["cost_per_m2"], m["carbon_per_m2"], brief.reuse_cycles)
    files = []
    widths = sorted(bom["modules"]) or [brief.reference_module_width]
    for width in widths:
        face = face_relief(best.screen, width, best.screen.height_above_desk)
        stem = f"module_{int(width)}x{int(best.screen.height_above_desk)}"
        write_svg(face, width, best.screen.height_above_desk, out / f"{stem}.svg", title=stem)
        write_obj(face, CORES[best.screen.core].thickness, out / f"{stem}.obj")
        files += [f"{stem}.svg", f"{stem}.obj"]
    write_bom_csv(bom, runs, out / "bill_of_materials.csv")
    with open(out / "design.json", "w", encoding="utf-8") as fh:
        json.dump({"best": design_json(best), "shortlist": [design_json(e) for e in shortlist],
                   "bill_of_materials": bom}, fh, indent=2, default=str)
    files += ["bill_of_materials.csv", "design.json"] + make_plots(best, history, out)
    write_report(out / "report.md", brief, best, shortlist, runs, bom, files, history)
    return bom


def cmd_run(args):
    brief = load_brief(args.brief)

    def progress(gen, best, mean):
        print(f"generation {gen:3d}  best {best:.3f}  mean {mean:.3f}")

    result = optimize.run(brief, args.population, args.generations, args.seed, progress=progress)
    shortlist = optimize.distinct(result.ranked)
    out = Path(args.out)
    bom = build_package(result.best, brief, out, result.history, shortlist)
    m = result.best.metrics
    print(f"\nBest: {result.best.screen}")
    print(f"alpha_w {m['alpha_w']:.2f} class {m['absorption_class']}, speech -{m['speech_reduction_dba']:.1f} dB(A), "
          f"heaviest module {m['heaviest_module_kg']:.1f} kg, {bom['module_count']} modules for the layout")
    print(f"Report: {out / 'report.md'}")


def cmd_evaluate(args):
    brief = load_brief(args.brief)
    with open(args.design, encoding="utf-8") as fh:
        data = json.load(fh)
    data = data.get("best", data).get("screen", data)
    e = evaluate(DeskScreen(**data), brief)
    build_package(e, brief, Path(args.out), [(0, e.fitness, e.fitness)], [e])
    print(json.dumps(e.metrics, indent=2, default=str))


def cmd_catalogue(_args):
    for title, table in (("Felts", FELTS), ("Cores", CORES), ("Joints", JOINTS), ("Supports", SUPPORTS)):
        print(title)
        for key, item in table.items():
            print(f"  {key:26s} {item.name}: {item.note}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="biothing_acoustics", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="breed a design for a brief")
    r.add_argument("--brief", default="briefs/exhibition_desk_screen.json")
    r.add_argument("--out", default="output")
    r.add_argument("--population", type=int, default=40)
    r.add_argument("--generations", type=int, default=30)
    r.add_argument("--seed", type=int, default=1)
    r.set_defaults(func=cmd_run)
    e = sub.add_parser("evaluate", help="evaluate a design JSON (e.g. an edited design.json)")
    e.add_argument("design")
    e.add_argument("--brief", default="briefs/exhibition_desk_screen.json")
    e.add_argument("--out", default="output/evaluated")
    e.set_defaults(func=cmd_evaluate)
    c = sub.add_parser("catalogue", help="list materials, joints and supports")
    c.set_defaults(func=cmd_catalogue)
    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
