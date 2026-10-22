"""Genome, decoding and evaluation of a desk screen.

A screen is a rigid core with a sculpted felt face on both sides. One felt
sheet is formed over a relief: where the relief is low the felt is pressed
(denser, higher flow resistivity, better at mid/high frequencies); where it
is high the felt keeps its thickness and stands off the core, leaving an air
cavity that pulls absorption down towards speech fundamentals. The cell
pattern and relief field spread these conditions over the face.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from . import acoustics, geometry, speech
from .materials import CORES, FELTS
from .modular import JOINTS, SUPPORTS, run_leak_fraction, tile

# (name, kind, choices or bounds). Genomes are vectors in [0, 1]^len(GENES).
GENES = (
    ("pattern", "choice", geometry.PATTERNS),
    ("cell_density", "float", (60.0, 400.0)),  # cells per m2
    ("field", "choice", geometry.FIELDS),
    ("field_exponent", "float", (0.4, 3.0)),
    ("felt", "choice", tuple(FELTS)),
    ("felt_thickness", "float", (6.0, 30.0)),  # mm, sheet before forming
    ("min_compression", "float", (0.3, 1.0)),  # thinnest formed felt / sheet thickness
    ("relief_depth", "float", (0.0, 60.0)),  # mm, deepest air cavity behind the felt
    ("core", "choice", tuple(CORES)),
    ("height_above_desk", "height", None),  # bounds from the brief
    ("joint", "choice", tuple(JOINTS)),
    ("support", "choice", tuple(SUPPORTS)),
)
STEPS = {"cell_density": 10, "felt_thickness": 0.5, "relief_depth": 0.5, "height_above_desk": 10,
         "field_exponent": 0.05, "min_compression": 0.01}


@dataclass(frozen=True)
class DeskScreen:
    pattern: str
    cell_density: float
    field: str
    field_exponent: float
    felt: str
    felt_thickness: float
    min_compression: float
    relief_depth: float
    core: str
    height_above_desk: float
    joint: str
    support: str


def decode(vector, brief) -> DeskScreen:
    values = {}
    for u, (name, kind, spec) in zip(vector, GENES):
        u = float(np.clip(u, 0.0, 1.0))
        if kind == "choice":
            allowed = getattr(brief, f"allowed_{name}s", None)
            if allowed:
                unknown = set(allowed) - set(spec)
                if unknown:
                    raise ValueError(f"unknown {name} keys in brief: {sorted(unknown)}")
                spec = tuple(k for k in spec if k in allowed)
            values[name] = spec[min(int(u * len(spec)), len(spec) - 1)]
            continue
        if kind == "height":
            spec = (brief.constraints.min_height_above_desk, brief.constraints.max_height_above_desk)
        lo, hi = spec
        step = STEPS[name]
        values[name] = round(round((lo + u * (hi - lo)) / step) * step, 4)
    return DeskScreen(**values)


def face_relief(screen, width, height):
    """Cells of one face plus per-cell relief height, felt thickness, gap and sigma (mm, Pa*s/m2)."""
    cells = geometry.build_cells(screen.pattern, float(width), float(height), int(screen.cell_density))
    felt = FELTS[screen.felt]
    t0 = screen.felt_thickness
    if felt.formable:
        u = geometry.relief_field(screen.field, cells, width, height) ** screen.field_exponent
        lowest, deepest = t0 * screen.min_compression, t0 + screen.relief_depth
    else:  # flat sheet
        u, lowest, deepest = np.ones(len(cells)), t0, t0
    relief = lowest + u * (deepest - lowest)
    thickness = np.minimum(t0, relief)
    gap = relief - thickness
    sigma = acoustics.compressed_sigma(felt.flow_resistivity, thickness / t0)
    area = np.array([c.area for c in cells])
    return {"cells": cells, "area": area, "relief": relief, "thickness": thickness, "gap": gap, "sigma": sigma}


@dataclass
class Evaluation:
    screen: DeskScreen
    metrics: dict
    curves: dict
    fitness: float

    def to_dict(self):
        return {"screen": asdict(self.screen), "fitness": self.fitness, "metrics": self.metrics}


def evaluate(screen: DeskScreen, brief) -> Evaluation:
    width, height = brief.reference_module_width, screen.height_above_desk
    face = face_relief(screen, width, height)
    absorption = acoustics.panel_absorption(face["area"] / face["area"].sum(), face["thickness"] / 1000,
                                            face["gap"] / 1000, face["sigma"])
    ap = acoustics.practical_alpha(absorption["random"])
    aw, shape, cls = acoustics.weighted_alpha([ap[f] for f in acoustics.RATING_OCTAVES])

    felt, core = FELTS[screen.felt], CORES[screen.core]
    felt_mass = felt.density * screen.felt_thickness / 1000  # kg/m2 per face, kept when pressed
    surface_mass = core.surface_mass + 2 * felt_mass
    thickness = core.thickness + 2 * float(face["relief"].max())
    heaviest = surface_mass * max(brief.module_widths) / 1000 * height / 1000

    joint, support = JOINTS[screen.joint], SUPPORTS[screen.support]
    run = tile(brief.scenario.reference_run_length, brief.module_widths)
    leak = run_leak_fraction(run, height, joint, support)
    tl_panel = acoustics.mass_law_tl(surface_mass, speech.OCTAVES)
    tl = acoustics.composite_tl(tl_panel, leak)
    top = brief.scenario.desk_height + height
    attenuation = speech.screen_attenuation(brief.scenario, top, sum(run), tl)
    priv = speech.privacy(brief.scenario, attenuation, brief.scenario.noise)

    cost_m2 = core.surface_mass * core.cost_per_kg + 2 * felt_mass * felt.cost_per_kg
    carbon_m2 = core.surface_mass * core.carbon_per_kg + 2 * felt_mass * felt.carbon_per_kg
    module_area = width * height / 1e6
    metrics = {
        "alpha_w": aw,
        "shape_indicator": shape,
        "absorption_class": cls,
        "nrc": acoustics.nrc(absorption["random"]),
        "practical_alpha": {str(k): v for k, v in ap.items()},
        "speech_band_alpha": float(np.mean([ap[500], ap[1000], ap[2000]])),
        "absorption_area_per_module_m2": 2 * module_area * acoustics.nrc(absorption["random"]),
        "cells_per_face": len(face["cells"]),
        "surface_mass_kg_m2": surface_mass,
        "thickness_mm": thickness,
        "reference_module_mass_kg": surface_mass * module_area,
        "heaviest_module_kg": heaviest,
        "leak_fraction": leak,
        "screen_top_above_floor_mm": top,
        "speech_reduction_dba": priv["speech_reduction_dba"],
        "speech_open_dba": priv["speech_open_dba"],
        "speech_screened_dba": priv["speech_screened_dba"],
        "sii_open": priv["sii_open"],
        "sii_screened": priv["sii_screened"],
        "privacy_gain": priv["privacy_gain"],
        "cost_per_m2": cost_m2,
        "carbon_per_m2": carbon_m2,
        "carbon_per_m2_per_event": carbon_m2 / max(1, brief.reuse_cycles),
        "felt_formed": felt.formable,
    }
    curves = {
        "freqs": absorption["freqs"], "alpha_random": absorption["random"], "alpha_normal": absorption["normal"],
        "speech_freqs": speech.OCTAVES, "attenuation": attenuation, "tl": tl, "tl_panel": tl_panel,
        "levels_open": priv["levels_open"], "levels_screened": priv["levels_screened"],
        "noise": brief.scenario.noise, "face": face,
    }
    return Evaluation(screen, metrics, curves, score(metrics, brief))


def score(m, brief) -> float:
    """Weighted sum of normalised objectives minus constraint penalties."""
    w, c = brief.weights, brief.constraints
    value = (w.absorption * m["alpha_w"]
             + w.speech_absorption * m["speech_band_alpha"]
             + w.privacy * m["privacy_gain"]
             + w.lightness * (1 - m["heaviest_module_kg"] / c.max_module_mass)
             + w.cost * (1 - m["cost_per_m2"] / brief.cost_reference_per_m2)
             + w.carbon * (1 - m["carbon_per_m2"] / brief.carbon_reference_per_m2))
    if m["heaviest_module_kg"] > c.max_module_mass:
        value -= 0.2 + (m["heaviest_module_kg"] / c.max_module_mass - 1)
    if m["thickness_mm"] > c.max_thickness:
        value -= 0.2 + (m["thickness_mm"] / c.max_thickness - 1)
    return float(value)
