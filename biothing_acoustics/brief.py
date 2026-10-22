"""Design brief: the scenario, constraints and priorities a screen is bred for."""
from __future__ import annotations

import json
from dataclasses import dataclass, field, fields

import numpy as np

from .speech import OCTAVES

# Typical busy exhibition hall, about 51 dB(A).
EXHIBITION_HALL_NOISE = {"125": 55, "250": 52, "500": 49, "1000": 46, "2000": 42, "4000": 38, "8000": 34}


@dataclass
class Scenario:
    desk_height: float = 740  # mm
    mouth_height: float = 1200  # seated talker, mm above floor
    ear_height: float = 1200
    talker_distance: float = 1000  # horizontal mouth-to-screen distance, mm
    listener_distance: float = 1000
    reference_run_length: float = 4200  # screen run used to score designs, mm
    background_noise: dict = field(default_factory=lambda: dict(EXHIBITION_HALL_NOISE))

    @property
    def noise(self):
        return np.array([float(self.background_noise[str(int(f))]) for f in OCTAVES])


@dataclass
class Constraints:
    max_module_mass: float = 8.0  # kg, one person carries and installs a module
    max_thickness: float = 120.0  # mm, including both sculpted faces
    min_height_above_desk: float = 350.0
    max_height_above_desk: float = 650.0


@dataclass
class Weights:
    absorption: float = 0.30  # weighted absorption coefficient alpha_w
    speech_absorption: float = 0.15  # mean practical alpha at 500-2000 Hz
    privacy: float = 0.30  # loss of speech intelligibility across the screen
    lightness: float = 0.10
    cost: float = 0.075
    carbon: float = 0.075


@dataclass
class Brief:
    name: str = "Exhibition desk screen"
    scenario: Scenario = field(default_factory=Scenario)
    constraints: Constraints = field(default_factory=Constraints)
    weights: Weights = field(default_factory=Weights)
    module_widths: list = field(default_factory=lambda: [600, 800, 1000, 1200, 1400, 1600])
    reference_module_width: float = 1200
    reuse_cycles: int = 10  # events each module is expected to serve
    cost_reference_per_m2: float = 60.0  # EUR/m2 of screen, materials only
    carbon_reference_per_m2: float = 15.0  # kgCO2e/m2 of screen
    layout: list = field(default_factory=list)
    # Optional whitelists of catalogue keys, e.g. PET felts and PET cores only.
    allowed_felts: list | None = None
    allowed_cores: list | None = None
    allowed_joints: list | None = None
    allowed_supports: list | None = None


def _build(cls, data):
    known = {f.name for f in fields(cls)}
    unknown = set(data) - known
    if unknown:
        raise ValueError(f"unknown {cls.__name__} keys: {sorted(unknown)}")
    return cls(**data)


def load_brief(path=None) -> Brief:
    if path is None:
        return Brief()
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    data.pop("_comment", None)
    for key, cls in (("scenario", Scenario), ("constraints", Constraints), ("weights", Weights)):
        if key in data:
            data[key] = _build(cls, data[key])
    return _build(Brief, data)
