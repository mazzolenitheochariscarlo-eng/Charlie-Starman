"""Fabrication and visualisation files: SVG forming pattern, OBJ mesh, CSV bill of materials."""
from __future__ import annotations

import csv

import numpy as np


def _shade(value):
    g = int(235 - 175 * value)
    return f"rgb({g},{g},{min(255, g + 12)})"


def write_svg(face, width, height, path, title=""):
    """Face pattern at 1:1 in mm. Darker cells stand further off the core (deeper cavity).

    Cell outlines are the forming lines for the thermoforming mould or CNC relief.
    """
    relief = face["relief"]
    span = np.ptp(relief) or 1.0
    caption = 14
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}mm" height="{height + caption}mm" '
           f'viewBox="0 0 {width} {height + caption}">',
           f'<rect x="0" y="0" width="{width}" height="{height}" fill="none" stroke="black" stroke-width="1"/>']
    for cell, h in zip(face["cells"], relief):
        pts = " ".join(f"{x:.2f},{height - y:.2f}" for x, y in cell.polygon)
        out.append(f'<polygon points="{pts}" fill="{_shade((h - relief.min()) / span)}" '
                   f'stroke="#333" stroke-width="0.4"><title>relief {h:.1f} mm</title></polygon>')
    out.append(f'<text x="2" y="{height + 10}" font-family="sans-serif" font-size="7">{title} - '
               f'relief {relief.min():.1f}-{relief.max():.1f} mm, {len(relief)} cells</text>')
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))


def write_obj(face, core_thickness, path):
    """Mesh of one module: core slab plus a prism per cell on both faces (z = thickness axis)."""
    verts, faces = [], []

    def prism(poly, z0, z1):
        base = len(verts) + 1
        k = len(poly)
        verts.extend((x, y, z0) for x, y in poly)
        verts.extend((x, y, z1) for x, y in poly)
        faces.append([base + i for i in range(k)][::-1])
        faces.append([base + k + i for i in range(k)])
        for i in range(k):
            j = (i + 1) % k
            faces.append([base + i, base + j, base + k + j, base + k + i])

    xs = [x for c in face["cells"] for x, _ in c.polygon]
    ys = [y for c in face["cells"] for _, y in c.polygon]
    rect = [(min(xs), min(ys)), (max(xs), min(ys)), (max(xs), max(ys)), (min(xs), max(ys))]
    half = core_thickness / 2
    prism(rect, -half, half)
    for cell, h in zip(face["cells"], face["relief"]):
        prism(cell.polygon, half, half + h)
        prism(cell.polygon, -half - h, -half)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# biothing desk screen module, units mm\n")
        fh.writelines(f"v {x:.2f} {y:.2f} {z:.2f}\n" for x, y, z in verts)
        fh.writelines("f " + " ".join(map(str, f)) + "\n" for f in faces)


def write_bom_csv(bom, runs, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["island", "run", "requested_mm", "modules_mm", "shortfall_mm"])
        for r in runs:
            w.writerow([r.island, r.kind, r.length, " + ".join(map(str, map(int, r.modules))), r.shortfall])
        w.writerow([])
        w.writerow(["module_width_mm", "count", "mass_each_kg"])
        for width, count in bom["modules"].items():
            w.writerow([width, count, bom["module_mass_kg"][width]])
        w.writerow([])
        w.writerow(["joints", bom["joints"]])
        w.writerow(["supports", bom["supports"]])
