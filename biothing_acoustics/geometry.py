"""Generative cell patterns for the sculpted felt face.

The patterns follow natural growth logics: phyllotaxis (sunflower heads),
relaxed Voronoi (foams, cell tissue) and hexagonal packing (honeycomb).
Each cell then gets a value from a relief field; design.py turns that value
into felt compression and air-gap depth, so every cell is tuned slightly
differently, like the multi-scale scales of a moth wing.

Units are millimetres, origin at the bottom-left corner of the panel face.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

GOLDEN_ANGLE = math.pi * (3 - math.sqrt(5))
PATTERNS = ("phyllotaxis", "voronoi", "hex")
FIELDS = ("radial", "vertical", "gyroid", "noise")
NEIGHBOURS = 24  # nearest seeds used to clip each Voronoi cell


@dataclass(frozen=True)
class Cell:
    polygon: tuple
    area: float
    centroid: tuple


def area_centroid(poly):
    a = cx = cy = 0.0
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        a += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    a *= 0.5
    if abs(a) < 1e-12:
        return 0.0, poly[0]
    return abs(a), (cx / (6 * a), cy / (6 * a))


def _clip(poly, mx, my, nx, ny):
    """Keep the part of a convex polygon on the seed side of a bisector."""
    out = []
    n = len(poly)
    for i in range(n):
        px, py = poly[i]
        qx, qy = poly[(i + 1) % n]
        dp = (px - mx) * nx + (py - my) * ny
        dq = (qx - mx) * nx + (qy - my) * ny
        if dp <= 0:
            out.append((px, py))
        if dp * dq < 0:
            t = dp / (dp - dq)
            out.append((px + t * (qx - px), py + t * (qy - py)))
    return out


def seed_points(pattern, width, height, density, seed=0):
    """Seed points for `density` cells per m2 on a width x height (mm) face."""
    n = max(4, int(round(density * width * height / 1e6)))
    if pattern == "phyllotaxis":
        c = math.sqrt(width * height / (math.pi * n))
        r_max = math.hypot(width, height) / 2
        pts, i = [], 0
        while c * math.sqrt(i) <= r_max:
            r, a = c * math.sqrt(i + 0.5), i * GOLDEN_ANGLE
            x, y = width / 2 + r * math.cos(a), height / 2 + r * math.sin(a)
            if 0 < x < width and 0 < y < height:
                pts.append((x, y))
            i += 1
        return np.array(pts)
    if pattern == "hex":
        s = math.sqrt(2 * width * height / (math.sqrt(3) * n))
        pts, row, y = [], 0, s * math.sqrt(3) / 4
        while y < height:
            x = s / 2 if row % 2 else s / 4
            while x < width:
                pts.append((x, y))
                x += s
            y += s * math.sqrt(3) / 2
            row += 1
        return np.array(pts)
    if pattern == "voronoi":
        rng = np.random.default_rng(seed)
        return rng.uniform((0, 0), (width, height), size=(n, 2))
    raise ValueError(f"unknown pattern {pattern!r}")


def voronoi_cells(points, width, height):
    pts = np.asarray(points, dtype=float)
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
    order = np.argsort(d2, axis=1)[:, 1:NEIGHBOURS + 1]
    rect = [(0.0, 0.0), (width, 0.0), (width, height), (0.0, height)]
    cells = []
    for i, (px, py) in enumerate(pts):
        poly = rect
        for j in order[i]:
            qx, qy = pts[j]
            poly = _clip(poly, (px + qx) / 2, (py + qy) / 2, qx - px, qy - py)
            if len(poly) < 3:
                break
        if len(poly) >= 3:
            area, centroid = area_centroid(poly)
            if area > 1e-6:
                cells.append(Cell(tuple(poly), area, centroid))
    return cells


@lru_cache(maxsize=2048)
def build_cells(pattern, width, height, density, seed=0):
    """Cells of a pattern; Voronoi seeds get two Lloyd relaxation steps (foam-like)."""
    cells = voronoi_cells(seed_points(pattern, width, height, density, seed), width, height)
    if pattern == "voronoi":
        for _ in range(2):
            cells = voronoi_cells([c.centroid for c in cells], width, height)
    return tuple(cells)


def relief_field(kind, cells, width, height, seed=0):
    """Value in [0, 1] for every cell; drives how deep the relief is there."""
    xy = np.array([c.centroid for c in cells])
    x, y = xy[:, 0], xy[:, 1]
    if kind == "radial":
        u = np.hypot(x - width / 2, y - height / 2) / (math.hypot(width, height) / 2)
    elif kind == "vertical":
        u = y / height
    elif kind == "gyroid":
        period, z = 400.0, math.pi / 4  # a planar slice through a gyroid surface
        gx, gy = 2 * np.pi * x / period, 2 * np.pi * y / period
        u = np.sin(gx) * np.cos(gy) + np.sin(gy) * math.cos(z) + math.sin(z) * np.cos(gx)
    elif kind == "noise":
        u = np.random.default_rng(seed + 7).uniform(size=len(cells))
    else:
        raise ValueError(f"unknown field {kind!r}")
    span = u.max() - u.min()
    return (u - u.min()) / span if span > 0 else np.full(len(cells), 0.5)
