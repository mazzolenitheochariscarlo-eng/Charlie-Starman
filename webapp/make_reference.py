"""Write reference results from the Python engine for webapp/engine.test.js."""
import json
from pathlib import Path

from biothing_acoustics import forms

CASES = {
    "spikes": {"depth": 120.0, "base": 0.3, "taper": 1.7, "porous": "pet_felt_medium"},
    "sponge": {"depth": 90.0, "outer": 0.3, "inner": 0.9, "gap": 0.25, "porous": "hemp_fibre"},
    "cloud": {"depth": 60.0, "outer": 0.4, "inner": 1.0, "porous": "sheep_wool_batt"},
    "bottles": {"depth": 140.0, "hole": 6.0, "plate": 8.0, "porosity": 0.04, "fill": 0.4, "ratio": 0.6,
                "porous": "cork_granulate", "rigid": "birch_plywood"},
    "coil": {"depth": 50.0, "width": 4.0, "cell": 50.0, "grading": 0.5, "coil": 0.9, "rigid": "pla_print"},
    "hybrid": {"depth": 150.0, "front": 0.5, "outer": 0.5, "inner": 1.0, "hole": 5.0, "plate": 4.0,
               "porosity": 0.08, "porous": "mycelium", "rigid": "rigid_pet"},
}

brief = forms.FormsBrief(target="exhibition")
out = {"weights": {t: list(forms.band_weights(t)) for t in forms.TARGETS}, "cases": {}}
for key, params in CASES.items():
    r = forms.assess(forms.FAMILIES[key], params, brief)
    out["cases"][key] = {"params": params, "alpha": [float(a) for a in r.alpha], "mass": r.mass, "score": r.score}
Path(__file__).with_name("reference.json").write_text(json.dumps(out, indent=1))
