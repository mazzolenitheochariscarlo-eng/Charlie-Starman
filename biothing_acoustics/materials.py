"""Material library for the sculpted felt faces and the screen core.

All numbers are indicative starting points for early-stage design studies.
Replace them with supplier data sheets (flow resistivity, density) and EPDs
(embodied carbon) before taking any product decision.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Felt:
    name: str
    flow_resistivity: float  # Pa*s/m2 at nominal (uncompressed) density
    density: float  # kg/m3 at nominal thickness
    cost_per_kg: float  # EUR/kg
    carbon_per_kg: float  # kgCO2e/kg, cradle-to-gate
    family: str  # recycling family, used for the mono-material check
    formable: bool  # can take a relief (thermoformed, moulded or grown)
    note: str = ""


@dataclass(frozen=True)
class Core:
    name: str
    thickness: float  # mm
    density: float  # kg/m3
    cost_per_kg: float
    carbon_per_kg: float
    family: str
    note: str = ""

    @property
    def surface_mass(self) -> float:
        """kg/m2"""
        return self.thickness / 1000 * self.density


FELTS = {
    "pet_felt_light": Felt("Recycled PET felt, light", 9_000, 35, 5.0, 1.6, "PET", True,
                           "Soft, very light; needs depth for low frequencies"),
    "pet_felt_medium": Felt("Recycled PET felt, medium", 20_000, 60, 5.0, 1.6, "PET", True,
                            "General-purpose thermoformable acoustic felt"),
    "pet_felt_dense": Felt("Recycled PET felt, dense", 45_000, 120, 5.5, 1.6, "PET", True,
                           "Robust surface, good at mid/high frequencies"),
    "wool_felt": Felt("Wool felt", 25_000, 100, 14.0, 10.0, "wool", False,
                      "Natural, excellent absorber; cut flat, not thermoformed"),
    "recycled_cotton": Felt("Recycled cotton fibre", 12_000, 45, 3.5, 0.8, "cotton", False,
                            "Low carbon; fire treatment usually required"),
    "mycelium": Felt("Mycelium composite", 60_000, 90, 9.0, 0.6, "biological", True,
                     "Grown in a mould: the most literal 'biothing'; check fire class and humidity"),
}

CORES = {
    "pet_board_9mm": Core("Rigid recycled PET board 9 mm", 9, 330, 6.0, 1.8, "PET",
                          "Mono-material with PET felt: fully recyclable together"),
    "pet_board_4mm": Core("High-density PET board 4 mm", 4, 900, 6.0, 1.8, "PET",
                          "Thin, heavier per volume, better sound reduction"),
    "birch_plywood_6mm": Core("Birch plywood 6 mm", 6, 680, 3.5, 0.7, "wood",
                              "Stiff, CNC-friendly, robust for repeated set-ups"),
    "honeycomb_cardboard_16mm": Core("Honeycomb cardboard 16 mm", 16, 50, 10.0, 1.1, "paper",
                                     "Ultra light for shipping; low sound reduction"),
    "acrylic_4mm": Core("Acrylic (PMMA) 4 mm", 4, 1190, 4.0, 6.0, "PMMA",
                        "Transparent core option; heavy and high carbon"),
}
