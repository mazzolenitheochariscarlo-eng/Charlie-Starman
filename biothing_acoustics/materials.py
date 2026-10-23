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
    "sheep_wool_batt": Felt("Sheep wool batt", 6_000, 25, 8.0, 1.0, "wool", False,
                            "Very soft and light; regulates humidity; needs a cover"),
    "hemp_fibre": Felt("Hemp fibre mat", 8_000, 40, 4.0, 0.5, "biological", False,
                       "Fast-growing plant, low carbon; slightly rough look"),
    "coconut_coir": Felt("Coconut coir mat", 3_000, 60, 3.0, 0.4, "biological", False,
                         "Coarse, open fibres: good front layer, weak alone"),
    "cork_granulate": Felt("Bonded cork granulate", 15_000, 110, 6.0, 0.3, "biological", True,
                           "Bark of the cork oak, harvested without felling; can be moulded"),
    "cellulose_fibre": Felt("Recycled paper cellulose", 20_000, 45, 1.5, 0.3, "paper", True,
                            "Very cheap and low carbon; can be moulded like papier-mache"),
    "wood_wool": Felt("Wood wool board", 1_500, 400, 1.5, 0.5, "wood", False,
                      "Wood shavings bound with cement or magnesite; robust, heavy, very open"),
    "mineral_wool": Felt("Stone wool", 15_000, 50, 2.5, 1.3, "mineral", False,
                         "Non-combustible reference absorber; must be fully covered"),
    "melamine_foam": Felt("Melamine foam", 10_000, 9, 15.0, 4.0, "plastic", True,
                          "Ultra light, cuttable into any shape; fossil-based"),
}

# Rigid materials for plates, resonator necks and channel walls.
RIGIDS = {
    "rigid_pet": Core("Rigid recycled PET", 0, 900, 6.0, 1.8, "PET", "Mono-material with PET felts"),
    "birch_plywood": Core("Birch plywood", 0, 680, 3.5, 0.7, "wood", "CNC-friendly, robust"),
    "solid_cardboard": Core("Solid recycled cardboard", 0, 700, 1.5, 0.9, "paper",
                            "Cheap, light, for short exhibitions"),
    "pla_print": Core("3D-printed PLA (plant starch)", 0, 1240, 25.0, 3.0, "PLA",
                      "Any shape, ideal for coiled channels; slow for large areas"),
    "fired_clay": Core("Fired clay", 0, 1900, 4.0, 0.4, "mineral",
                       "Like the resonating pots built into ancient theatres and medieval churches"),
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
