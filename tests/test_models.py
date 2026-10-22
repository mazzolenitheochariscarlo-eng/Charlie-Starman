import json

import numpy as np
import pytest

from biothing_acoustics import acoustics, geometry, optimize, speech
from biothing_acoustics.brief import load_brief
from biothing_acoustics.cli import main
from biothing_acoustics.design import GENES, decode, evaluate
from biothing_acoustics.modular import JOINTS, SUPPORTS, plan_layout, run_leak_fraction, tile


def absorption(thickness, gap=0.0, sigma=20_000):
    return acoustics.panel_absorption([1.0], [thickness], [gap], [sigma])["random"]


def test_absorption_is_bounded_and_physical():
    thin, thick = absorption(0.01), absorption(0.04)
    assert np.all((thin >= 0) & (thin <= 1))
    i500 = list(acoustics.THIRD_OCTAVES).index(500)
    assert thick[i500] > thin[i500]
    with_gap = absorption(0.02, gap=0.05)
    assert with_gap[i500] > absorption(0.02)[i500]


def test_compression_raises_flow_resistivity():
    assert acoustics.compressed_sigma(20_000, 0.5) == pytest.approx(20_000 * 0.5 ** -1.53)


def test_weighted_alpha_iso_11654():
    assert acoustics.weighted_alpha([0.30, 0.60, 0.80, 0.90, 0.90]) == (0.60, "H", "C")
    assert acoustics.weighted_alpha([1, 1, 1, 1, 1]) == (1.0, "", "A")
    assert acoustics.weighted_alpha([0.0] * 5)[2] == "not classified"


def test_mass_law_and_leaks():
    assert acoustics.mass_law_tl(10, 500) == pytest.approx(20 * np.log10(5000) - 47)
    assert acoustics.composite_tl(40.0, 0.01) == pytest.approx(-10 * np.log10(0.99e-4 + 0.01))


def test_barrier_attenuation():
    att = speech.barrier_attenuation(np.array([-0.5, 0.0, 1.0, 1000.0]))
    assert att[0] == 0
    assert att[1] == pytest.approx(5.0)
    assert 12.5 < att[2] < 14
    assert att[3] == speech.MAX_BARRIER_ATTENUATION
    assert speech.BAND_IMPORTANCE.sum() == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("pattern", geometry.PATTERNS)
def test_cells_tile_the_panel(pattern):
    cells = geometry.build_cells(pattern, 1200.0, 500.0, 200)
    assert sum(c.area for c in cells) == pytest.approx(1200 * 500, rel=1e-6)
    assert 0.75 * 120 < len(cells) < 1.25 * 120


def test_tiling_and_layout():
    widths = [600, 800, 1000, 1200, 1400, 1600]
    assert tile(4200, widths) == (1400, 1400, 1400)
    assert tile(800, widths) == (800,)
    assert sum(tile(700, widths)) == 600
    runs = plan_layout([{"type": "bench", "name": "A", "desks_per_side": 3, "desk_width": 1400,
                         "desk_depth": 800}], widths)
    assert [r.kind for r in runs].count("divider") == 4
    leak = run_leak_fraction((1200, 1200), 500, JOINTS["slot_key"], SUPPORTS["desk_saddle"])
    assert leak == pytest.approx((3 * 500 + 10 * 2400) / (2400 * 510))


def test_decode_respects_brief_whitelists():
    brief = load_brief("briefs/pet_mono_material.json")
    rng = np.random.default_rng(0)
    for _ in range(50):
        s = decode(rng.uniform(size=len(GENES)), brief)
        assert s.felt.startswith("pet_") and s.core.startswith("pet_")
        assert 350 <= s.height_above_desk <= 650


def test_evaluate_and_optimise():
    brief = load_brief("briefs/exhibition_desk_screen.json")
    e = evaluate(decode(np.full(len(GENES), 0.5), brief), brief)
    assert 0 <= e.metrics["alpha_w"] <= 1
    assert e.metrics["speech_screened_dba"] < e.metrics["speech_open_dba"]
    result = optimize.run(brief, population=8, generations=2, seed=3)
    assert result.history[-1][1] >= result.history[0][1]


def test_cli_writes_package(tmp_path):
    main(["run", "--brief", "briefs/exhibition_desk_screen.json", "--out", str(tmp_path),
          "--population", "6", "--generations", "1"])
    assert (tmp_path / "report.md").exists()
    assert list(tmp_path.glob("module_*.svg")) and list(tmp_path.glob("module_*.obj"))
    data = json.loads((tmp_path / "design.json").read_text())
    main(["evaluate", str(tmp_path / "design.json"), "--out", str(tmp_path / "again")])
    again = json.loads((tmp_path / "again" / "design.json").read_text())
    assert again["best"]["metrics"]["alpha_w"] == data["best"]["metrics"]["alpha_w"]
