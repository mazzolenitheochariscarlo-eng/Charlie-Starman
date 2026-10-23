import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from biothing_acoustics import forms, layers

ROOT = Path(__file__).resolve().parents[1]


def test_band_weights_sum_to_one():
    for target in forms.TARGETS:
        w = forms.band_weights(target)
        assert w.sum() == pytest.approx(1.0)
        assert len(w) == len(forms.FREQS)


def test_limp_frame_limits():
    f = np.array([[250.0]])
    rigid, _ = layers.porous_medium(f, 10_000)
    heavy, _ = layers.porous_medium(f, 10_000, frame_density=1e9)
    weightless, _ = layers.porous_medium(f, 10_000, frame_density=0.0)
    assert heavy[0, 0] == pytest.approx(rigid[0, 0], rel=1e-6)
    assert weightless[0, 0] == pytest.approx(layers.RHO0, rel=1e-6)


def test_free_layer_conserves_energy():
    f = forms.FREQS
    m = layers.layer(f, forms.THETA, *layers.porous_medium(f, 20_000), 0.05)
    assert np.all((layers.alpha_free(m, forms.THETA) >= 0) & (layers.alpha_free(m, forms.THETA) <= 1))


@pytest.mark.parametrize("key", list(forms.FAMILIES))
def test_every_family_breeds(key):
    brief = forms.FormsBrief(max_depth=100, mass_budget=3)
    best, history = forms.breed(forms.FAMILIES[key], brief, population=8, generations=3)
    limit = 2 if forms.FAMILIES[key].mounting == "hanging" else 1
    assert np.all((best.alpha >= 0) & (best.alpha <= limit))
    assert history[-1][1] >= history[0][1]


def test_deeper_absorbs_more_at_low_pitch():
    shallow = forms.assess(forms.FAMILIES["sponge"], {"depth": 30, "outer": 1, "inner": 1, "gap": 0, "porous": "pet_felt_medium"},
                           forms.FormsBrief())
    deep = forms.assess(forms.FAMILIES["sponge"], {"depth": 150, "outer": 1, "inner": 1, "gap": 0, "porous": "pet_felt_medium"},
                        forms.FormsBrief())
    assert deep.alpha[4] > shallow.alpha[4]  # 250 Hz


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_web_engines_match_python():
    subprocess.run(["python3", "webapp/make_reference.py"], cwd=ROOT, check=True, env={"PYTHONPATH": str(ROOT), "PATH": "/usr/bin:/bin"})
    for test in ("webapp/engine.test.js", "webapp/rooms.test.js", "webapp/growth.test.js"):
        subprocess.run(["node", test], cwd=ROOT, check=True)
