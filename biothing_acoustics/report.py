"""Markdown report and plots for a pipeline run."""
from __future__ import annotations

from dataclasses import asdict

from .materials import CORES, FELTS
from .modular import JOINTS, SUPPORTS


def make_plots(evaluation, history, out_dir):
    """PNG plots; returns file names, or an empty list if matplotlib is missing."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
    except ImportError:
        return []
    c, files = evaluation.curves, []

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogx(c["freqs"], c["alpha_random"], "o-", label="random incidence (predicted)")
    ax.semilogx(c["freqs"], c["alpha_normal"], "--", label="normal incidence")
    ap = evaluation.metrics["practical_alpha"]
    ax.step([int(k) for k in ap], list(ap.values()), where="mid", label="practical alpha_p (octave)")
    ax.set(xlabel="Frequency (Hz)", ylabel="Absorption coefficient", ylim=(0, 1.05),
           title=f"alpha_w = {evaluation.metrics['alpha_w']:.2f} (class {evaluation.metrics['absorption_class']})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "absorption.png", dpi=150)
    files.append("absorption.png")

    fig, ax = plt.subplots(figsize=(7, 4))
    f = c["speech_freqs"]
    ax.semilogx(f, c["levels_open"], "o-", label="speech, no screen")
    ax.semilogx(f, c["levels_screened"], "o-", label="speech, with screen")
    ax.semilogx(f, c["noise"], "k--", label="background noise")
    ax.set(xlabel="Octave band (Hz)", ylabel="dB", title="Speech at the facing desk")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "speech_privacy.png", dpi=150)
    files.append("speech_privacy.png")

    fig, ax = plt.subplots(figsize=(7, 3.5))
    gens, best, mean = zip(*history)
    ax.plot(gens, best, label="best")
    ax.plot(gens, mean, label="population mean")
    ax.set(xlabel="Generation", ylabel="Fitness", title="Evolution")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "evolution.png", dpi=150)
    files.append("evolution.png")

    face = c["face"]
    fig, ax = plt.subplots(figsize=(8, 8 * evaluation.screen.height_above_desk / 1200 + 0.8))
    coll = PolyCollection([cell.polygon for cell in face["cells"]], array=face["relief"],
                          cmap="viridis", edgecolors="white", linewidths=0.4)
    ax.add_collection(coll)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set(xlabel="mm", ylabel="mm", title="Face relief (mm from core)")
    fig.colorbar(coll, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(out_dir / "face_pattern.png", dpi=150)
    files.append("face_pattern.png")
    plt.close("all")
    return files


def _design_table(e):
    s = e.screen
    if e.metrics["felt_formed"]:
        rows = [
            ("Cell pattern", f"{s.pattern}, {s.cell_density:.0f} cells/m2 ({e.metrics['cells_per_face']} per face)"),
            ("Relief field", f"{s.field}, exponent {s.field_exponent}"),
            ("Felt", f"{FELTS[s.felt].name}, sheet {s.felt_thickness} mm"),
            ("Forming", f"pressed to {s.min_compression:.0%} of sheet thickness, cavities up to {s.relief_depth} mm"),
        ]
    else:
        rows = [("Felt", f"{FELTS[s.felt].name}, flat sheet {s.felt_thickness} mm (not formable, no relief)")]
    rows += [
        ("Core", CORES[s.core].name),
        ("Height above desk", f"{s.height_above_desk:.0f} mm (top at {e.metrics['screen_top_above_floor_mm']:.0f} mm)"),
        ("Total thickness", f"{e.metrics['thickness_mm']:.0f} mm"),
        ("Joint", JOINTS[s.joint].name),
        ("Support", SUPPORTS[s.support].name),
    ]
    return "\n".join(f"| {k} | {v} |" for k, v in rows)


def write_report(path, brief, best, shortlist, runs, bom, files, history):
    m = best.metrics
    ap = m["practical_alpha"]
    lines = [
        f"# {brief.name}: generated design",
        "",
        f"Bred over {len(history) - 1} generations; best fitness {best.fitness:.3f}.",
        "",
        "## Winning design",
        "",
        "| Parameter | Value |",
        "|---|---|",
        _design_table(best),
        "",
        "## Acoustic performance",
        "",
        f"- Weighted absorption **alpha_w = {m['alpha_w']:.2f}{'(' + m['shape_indicator'] + ')' if m['shape_indicator'] else ''}**, "
        f"class **{m['absorption_class']}**, NRC {m['nrc']:.2f}",
        f"- Equivalent absorption area per {brief.reference_module_width:.0f} mm module (both faces): "
        f"{m['absorption_area_per_module_m2']:.2f} m2 Sabine",
        f"- Speech at the facing desk: {m['speech_open_dba']:.1f} dB(A) without screen, "
        f"{m['speech_screened_dba']:.1f} dB(A) with screen (-{m['speech_reduction_dba']:.1f} dB)",
        f"- Speech intelligibility index: {m['sii_open']:.2f} -> {m['sii_screened']:.2f} "
        f"(privacy gain {m['privacy_gain']:.0%})",
        f"- Open area from joints and the slit above the desk: {m['leak_fraction']:.1%}",
        "",
        "| Octave (Hz) | " + " | ".join(ap) + " |",
        "|---" * (len(ap) + 1) + "|",
        "| alpha_p | " + " | ".join(f"{v:.2f}" for v in ap.values()) + " |",
        "",
        "## Exhibition system",
        "",
        "| Item | Value |",
        "|---|---|",
        "| Modules | " + ", ".join(f"{c} x {w} mm" for w, c in bom["modules"].items()) + " |",
        "| Module mass | " + ", ".join(f"{w} mm: {kg} kg" for w, kg in bom["module_mass_kg"].items()) + " |",
        f"| Joints / supports | {bom['joints']} / {bom['supports']} |",
        f"| Screen area | {bom['screen_area_m2']:.1f} m2 |",
        f"| Total mass incl. fittings | {bom['total_mass_kg']:.0f} kg |",
        f"| Materials and fittings cost | {bom['cost_eur']:.0f} EUR |",
        f"| Embodied carbon | {bom['carbon_kg']:.0f} kgCO2e ({bom['carbon_per_event_kg']:.0f} per event over {brief.reuse_cycles} events) |",
        f"| Flat-pack volume | {bom['flatpack_volume_m3']:.2f} m3 (~{bom['pallets']} pallet(s)) |",
        f"| Assembly effort | {bom['assembly_person_hours']:.1f} person-hours |",
        f"| Mono-material (felt and core recycle together) | {'yes' if bom['mono_material'] else 'no'} |",
        "",
    ]
    if bom["shortfalls_mm"]:
        lines += ["Runs shorter than requested (normal for dividers; check spines and fronts):", ""]
        lines += [f"- {island} {kind}: {gap:.0f} mm short" for island, kind, gap in bom["shortfalls_mm"]]
        lines.append("")
    lines += ["## Shortlist (distinct alternatives)", "",
              "| # | Pattern | Felt | Core | Height | alpha_w | Speech -dB(A) | Heaviest module | Fitness |",
              "|---|---|---|---|---|---|---|---|---|"]
    for i, e in enumerate(shortlist, 1):
        s, mm = e.screen, e.metrics
        pattern = s.pattern if mm["felt_formed"] else "flat"
        lines.append(f"| {i} | {pattern} | {s.felt} | {s.core} | {s.height_above_desk:.0f} | "
                     f"{mm['alpha_w']:.2f} | {mm['speech_reduction_dba']:.1f} | {mm['heaviest_module_kg']:.1f} kg | "
                     f"{e.fitness:.3f} |")
    lines += ["", "## Files", ""] + [f"- `{f}`" for f in files]
    lines += ["", "## Limits of the models", "",
              "- Absorption is predicted for an infinite surface (Miki model, Paris integration). Measured "
              "values on real panels (ISO 354) are usually higher because of edge effects.",
              "- Speech privacy counts the direct, diffracted and transmitted paths only. Ceiling and wall "
              "reflections are not modelled.",
              "- Material, cost and carbon values are indicative. Replace them with supplier data and EPDs.",
              "- Confirm with impedance-tube samples (ISO 10534-2) and then a reverberation-room test.", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def design_json(evaluation):
    return {"screen": asdict(evaluation.screen), "fitness": evaluation.fitness, "metrics": evaluation.metrics}
