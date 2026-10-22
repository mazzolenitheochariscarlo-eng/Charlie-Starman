"""Modular screen system for temporary exhibitions.

Design rules behind the catalogue:
- dry assembly only (no glue, no screws into the panels), so modules can be
  set up, taken down and reused for many events;
- a small family of module widths that tile common desk sizes;
- every module light enough for one person to carry and install;
- flat-pack transport.
Joint and support data are indicative; replace them with real hardware specs.
"""
from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

from .materials import CORES, FELTS


@dataclass(frozen=True)
class Joint:
    name: str
    gap_mm: float  # open slit left between two modules
    mass_kg: float
    cost: float  # EUR each
    carbon_kg: float
    assembly_s: float
    note: str


@dataclass(frozen=True)
class Support:
    name: str
    bottom_gap_mm: float  # open slit between desk top and screen
    mass_kg: float
    cost: float
    carbon_kg: float
    assembly_s: float
    note: str


JOINTS = {
    "slot_key": Joint("Slot and key, CNC-cut", 3, 0.05, 1.2, 0.15, 20,
                      "No tools, no glue; key cut from the core offcuts"),
    "magnetic_strip": Joint("Magnetic edge strip", 1, 0.12, 4.5, 0.6, 10,
                            "Fastest and best seal; separate magnets at end of life"),
    "hinge_clip": Joint("Hinge clip", 4, 0.09, 2.8, 0.4, 25,
                        "Any angle: zigzag or curved runs, can stand without clamps"),
}

SUPPORTS = {
    "desk_clamp": Support("Desk-edge clamp", 20, 0.45, 8.0, 2.5, 60,
                          "Most stable; leaves a gap above the desk"),
    "weighted_foot": Support("Weighted foot", 5, 1.2, 6.0, 2.0, 15,
                             "Free-standing on the desk, no tools; heavier to ship"),
    "desk_saddle": Support("Folded board saddle", 10, 0.25, 3.5, 0.5, 20,
                           "Straddles the joint between facing desks; same material as the core"),
}

MODULE_HANDLING_S = 30  # unpack, position and level one module
PALLET_VOLUME_M3 = 1.2 * 0.8 * 1.6 * 0.7  # euro pallet, 1.6 m high, 70 % fill


@dataclass(frozen=True)
class Run:
    island: str
    kind: str  # spine, front or divider
    length: float  # mm requested
    modules: tuple  # module widths in mm

    @property
    def shortfall(self) -> float:
        return self.length - sum(self.modules)


def tile(length, widths, step=50):
    """Fewest modules from `widths` that fill `length` (mm); ties favour equal widths.

    If no exact fit exists, the longest fill shorter than `length` is returned.
    """
    target = int(round(length / step))
    sizes = sorted({int(round(w / step)) for w in widths})
    best = [None] * (target + 1)
    best[0] = ()
    for x in range(1, target + 1):
        for w in sizes:
            prev = best[x - w] if w <= x else None
            if prev is None:
                continue
            cand = tuple(sorted(prev + (w,), reverse=True))
            cur = best[x]
            if cur is None or (len(cand), cand[0] - cand[-1]) < (len(cur), cur[0] - cur[-1]):
                best[x] = cand
    for x in range(target, -1, -1):
        if best[x] is not None:
            return tuple(w * step for w in best[x])
    return ()


def run_leak_fraction(modules, height, joint, support):
    """Open area fraction of a run: joints between modules plus the slit under it."""
    length = sum(modules)
    if length == 0:
        return 0.0
    leak = (len(modules) - 1) * joint.gap_mm * height + support.bottom_gap_mm * length
    return leak / (length * (height + support.bottom_gap_mm))


def plan_layout(layout, widths):
    """Screen runs for a list of desk islands.

    bench: desks facing each other -> one central spine plus dividers between
           neighbours on both sides.
    row:   a single row of desks facing visitors -> one front run plus dividers.
    """
    runs = []
    for item in layout:
        n = item["desks_per_side"] if item["type"] == "bench" else item["desks"]
        dw, dd = item["desk_width"], item["desk_depth"]
        count = item.get("count", 1)
        for k in range(count):
            name = item["name"] + (f" #{k + 1}" if count > 1 else "")
            main = "spine" if item["type"] == "bench" else "front"
            runs.append(Run(name, main, n * dw, tile(n * dw, widths)))
            if item.get("dividers", True):
                dividers = 2 * (n - 1) if item["type"] == "bench" else n - 1
                for _ in range(dividers):
                    runs.append(Run(name, "divider", dd, tile(dd, widths)))
    return runs


def bill_of_materials(screen, runs, surface_mass, thickness, cost_m2, carbon_m2, reuse_cycles):
    """Counts, mass, cost, carbon, transport and assembly effort for a layout."""
    joint, support = JOINTS[screen.joint], SUPPORTS[screen.support]
    h = screen.height_above_desk
    modules = Counter(w for run in runs for w in run.modules)
    n_modules = sum(modules.values())
    area = sum(w * h for w, c in modules.items() for _ in range(c)) / 1e6
    n_joints = sum(max(0, len(r.modules) - 1) for r in runs) + sum(r.kind == "divider" for r in runs)
    n_supports = sum(len(r.modules) + 1 for r in runs if r.modules)
    panel_mass = area * surface_mass
    fittings_mass = n_joints * joint.mass_kg + n_supports * support.mass_kg
    cost = area * cost_m2 + n_joints * joint.cost + n_supports * support.cost
    carbon = area * carbon_m2 + n_joints * joint.carbon_kg + n_supports * support.carbon_kg
    volume = area * thickness / 1000 * 1.15
    assembly_s = n_modules * MODULE_HANDLING_S + n_joints * joint.assembly_s + n_supports * support.assembly_s
    return {
        "modules": dict(sorted(modules.items())),
        "module_count": n_modules,
        "module_mass_kg": {w: round(w * h / 1e6 * surface_mass, 2) for w in sorted(modules)},
        "joints": n_joints,
        "supports": n_supports,
        "screen_area_m2": area,
        "total_mass_kg": panel_mass + fittings_mass,
        "cost_eur": cost,
        "carbon_kg": carbon,
        "carbon_per_event_kg": carbon / max(1, reuse_cycles),
        "flatpack_volume_m3": volume,
        "pallets": max(1, math.ceil(volume / PALLET_VOLUME_M3)),
        "assembly_person_hours": assembly_s / 3600,
        "shortfalls_mm": [(r.island, r.kind, r.shortfall) for r in runs if r.shortfall > 0],
        "mono_material": FELTS[screen.felt].family == CORES[screen.core].family,
    }
