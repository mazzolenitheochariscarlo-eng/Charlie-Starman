"""Forms lab: breed and compare general absorber shapes, independent of any product.

Six shape families, each imitating a natural strategy, are given the same
budget (maximum depth, material mass per m2) and the same target sound, then
each is evolved separately so the comparison is best-against-best.

Scores are absorption per m2 of footprint. Wall forms are backed by a hard
wall. The hanging cloud is open on both sides, so it is scored per m2 of
hanging element for both faces and can exceed 1.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import layers as L
from .acoustics import RHO0, THIRD_OCTAVES
from .materials import FELTS, RIGIDS
from .optimize import evolve

THETA = np.linspace(0, np.radians(89), 46)[None, :]
FREQS = THIRD_OCTAVES[:, None]
N_SUB = 14  # slices used to describe graded shapes
WALL = 0.002  # m, thickness of channel walls and partitions
PARTITION_SPACING = 0.1  # m, cavity partitions behind resonator plates

# How much each octave matters for a target (weights sum to 1).
TARGETS = {
    "voices": {125: 0.02, 250: 0.10, 500: 0.20, 1000: 0.28, 2000: 0.25, 4000: 0.15},
    "crowd": {125: 0.15, 250: 0.25, 500: 0.25, 1000: 0.20, 2000: 0.10, 4000: 0.05},
    "exhibition": {125: 0.085, 250: 0.175, 500: 0.225, 1000: 0.24, 2000: 0.175, 4000: 0.10},
    "low_noise": {125: 0.35, 250: 0.30, 500: 0.20, 1000: 0.10, 2000: 0.05, 4000: 0.0},
}


def band_weights(target):
    octaves = np.array(sorted(TARGETS[target]))
    nearest = octaves[np.argmin(np.abs(np.log2(THIRD_OCTAVES[:, None] / octaves[None, :])), axis=1)]
    w = np.array([TARGETS[target][o] for o in nearest]) / np.bincount(
        np.searchsorted(octaves, nearest))[np.searchsorted(octaves, nearest)]
    return w / w.sum()


@dataclass
class FormsBrief:
    max_depth: float = 200.0  # mm, room you allow the absorber to take
    mass_budget: float = 6.0  # kg per m2 of footprint
    target: str = "exhibition"
    porous: list | None = None  # allowed FELTS keys (None = all)
    rigid: list | None = None  # allowed RIGIDS keys (None = all)
    families: list | None = None


# ---------------------------------------------------------------- models ----

def _rigid_admittance_alpha(y):
    return L.paris(L.alpha_from_admittance(y, THETA), THETA[0])


def graded_matrix(sigma, fractions, depth):
    dz = depth / len(fractions)
    return L.stack([L.layer(FREQS, THETA, *L.mixed_medium(FREQS, sigma, s), dz) for s in fractions])


def air_matrix(depth, theta=THETA):
    return L.layer(FREQS, theta, RHO0 + 0j, L.K0 + 0j, depth)


def _mid():
    return (np.arange(N_SUB) + 0.5) / N_SUB


def spikes(p, brief):
    """Wedges / spikes: solid base, then a taper that thins to a point."""
    felt = FELTS[p["porous"]]
    d, base = p["depth"] / 1000, p["base"]
    z = _mid()  # 0 = tip, 1 = wall
    taper = 1 - base
    s = np.where(z >= taper, 1.0, (z / taper) ** p["taper"]) if taper > 0 else np.ones(N_SUB)
    alpha = _rigid_admittance_alpha(L.admittance_on(graded_matrix(felt.flow_resistivity, s, d), 0))
    return alpha, s.mean() * d * felt.density, {"fractions": s}


def sponge(p, brief):
    """Gradient sponge: open at the surface, denser towards the wall, optional air gap behind."""
    felt = FELTS[p["porous"]]
    d = p["depth"] / 1000
    d_sponge = d * (1 - p["gap"])
    s = p["outer"] + (p["inner"] - p["outer"]) * _mid()
    m = graded_matrix(felt.flow_resistivity, s, d_sponge)
    if p["gap"] > 0:
        m = L.multiply(m, air_matrix(d * p["gap"]))
    alpha = _rigid_admittance_alpha(L.admittance_on(m, 0))
    return alpha, s.mean() * d_sponge * felt.density, {"fractions": s}


def cloud(p, brief):
    """Hanging cloud: a sponge hung in the air, open on both faces (scored for both)."""
    felt = FELTS[p["porous"]]
    d = p["depth"] / 1000
    half = p["outer"] + (p["inner"] - p["outer"]) * _mid()
    s = np.concatenate([half, half[::-1]])
    dz = d / len(s)
    m = L.stack([L.layer(FREQS, THETA, *L.mixed_medium(FREQS, felt.flow_resistivity, x, felt.density), dz)
                 for x in s])
    alpha = 2 * L.paris(L.alpha_free(m, THETA), THETA[0])
    return alpha, s.mean() * d * felt.density, {"fractions": s}


def _resonator_impedance(f, p, depth, fill, sigma):
    """Perforated plate in front of a cavity whose front part is filled with felt."""
    plate = L.perforated_plate(f, p["hole"] / 1000, p["plate"] / 1000, p["porosity"])
    felt = L.layer(f, 0.0, *L.porous_medium(f, sigma), depth * fill)
    air = L.layer(f, 0.0, RHO0 + 0j, L.K0 + 0j, depth * (1 - fill))
    y = L.admittance_on(L.multiply(felt, air), 0)
    return plate + 1 / y


def bottles(p, brief):
    """Bottle resonators: holes in a plate, cavities behind, two cavity sizes side by side."""
    felt, rigid = FELTS[p["porous"]], RIGIDS[p["rigid"]]
    d = p["depth"] / 1000
    depths = (d, d * p["ratio"])
    y = sum(0.5 / _resonator_impedance(FREQS, p, x, p["fill"], felt.flow_resistivity) for x in depths)
    alpha = _rigid_admittance_alpha(y)  # partitioned cavities react locally
    mean_depth = np.mean(depths)
    mass = (p["plate"] / 1000 * (1 - p["porosity"]) * rigid.density
            + p["fill"] * mean_depth * felt.density
            + 2 / PARTITION_SPACING * mean_depth * WALL * rigid.density)
    return alpha, mass, {"depths": depths}


def coil(p, brief):
    """Coiled channels: long narrow tubes folded into a thin slab, four lengths like a cochlea."""
    rigid = RIGIDS[p["rigid"]]
    d, w, a = p["depth"] / 1000, p["width"] / 1000, p["cell"] / 1000
    l_max = p["coil"] * a * d / (w + WALL)
    lengths = l_max * (1 - p["grading"] * np.arange(4) / 3)
    eps = w / a
    rho, k = L.slit_medium(FREQS, w)
    zt, kt = np.sqrt(rho * k), 2 * np.pi * FREQS * np.sqrt(rho / k)
    end = 1j * 2 * np.pi * FREQS * RHO0 * 0.6 * w
    y = sum(0.25 * eps / (-1j * zt / np.tan(kt * length) + end) for length in lengths)
    alpha = _rigid_admittance_alpha(y)
    mass = (WALL / (w + WALL) * d * p["coil"] + 0.002 * (1 - eps)) * rigid.density
    return alpha, mass, {"lengths": lengths}


def hybrid(p, brief):
    """Moss on a shell: a gradient sponge in front of a bottle-resonator layer."""
    felt, rigid = FELTS[p["porous"]], RIGIDS[p["rigid"]]
    d = p["depth"] / 1000
    front = d * p["front"]
    s = p["outer"] + (p["inner"] - p["outer"]) * _mid()
    cavity = max(d - front - p["plate"] / 1000, 0.001)
    z_res = _resonator_impedance(FREQS, p, cavity, 0.0, felt.flow_resistivity)
    alpha = _rigid_admittance_alpha(L.admittance_on(graded_matrix(felt.flow_resistivity, s, front), 1 / z_res))
    mass = (s.mean() * front * felt.density + p["plate"] / 1000 * (1 - p["porosity"]) * rigid.density
            + 2 / PARTITION_SPACING * cavity * WALL * rigid.density)
    return alpha, mass, {"fractions": s, "cavity": cavity}


# ------------------------------------------------------------- families ----

@dataclass(frozen=True)
class Family:
    key: str
    name: str
    imitates: str
    idea: str
    genes: tuple  # (name, kind, spec); kinds: float, depth, porous, rigid
    model: object
    mounting: str = "wall"


DEPTH = ("depth", "depth", (0.15, 1.0))
FAMILIES = {f.key: f for f in (
    Family("spikes", "Spikes and wedges", "grass, fur, owl-feather velvet",
           "Sound slides in between the points without a hard surface to bounce off.",
           (DEPTH, ("base", "float", (0.0, 0.7)), ("taper", "float", (0.4, 3.0)), ("porous", "porous", None)),
           spikes),
    Family("sponge", "Gradient sponge", "sea sponge, moss, bone",
           "Open and airy at the surface so sound gets in, denser inside so it gets trapped.",
           (DEPTH, ("outer", "float", (0.05, 1.0)), ("inner", "float", (0.05, 1.0)), ("gap", "float", (0.0, 0.6)),
            ("porous", "porous", None)),
           sponge),
    Family("cloud", "Hanging cloud", "hanging moss, clouds of seeds, a forest canopy",
           "Hangs free in the air, so sound is caught from both sides.",
           (DEPTH, ("outer", "float", (0.05, 1.0)), ("inner", "float", (0.05, 1.0)), ("porous", "porous", None)),
           cloud, "hanging"),
    Family("bottles", "Bottle resonators", "seashells, gourds, hollow stems",
           "Air in a small neck bounces on the air in the cavity and swallows one low pitch.",
           (DEPTH, ("hole", "float", (2.0, 20.0)), ("plate", "float", (2.0, 15.0)), ("porosity", "float", (0.005, 0.2)),
            ("fill", "float", (0.0, 1.0)), ("ratio", "float", (0.2, 1.0)), ("porous", "porous", None),
            ("rigid", "rigid", None)),
           bottles),
    Family("coil", "Coiled channels", "snail shell, the cochlea in your ear",
           "A very long, narrow tube folded into a thin object catches low pitches without needing depth.",
           (DEPTH, ("width", "float", (2.0, 15.0)), ("cell", "float", (20.0, 80.0)), ("grading", "float", (0.0, 0.9)),
            ("coil", "float", (0.3, 1.0)), ("rigid", "rigid", None)),
           coil),
    Family("hybrid", "Moss on a shell", "moss growing on a shell or bark",
           "A soft sponge for high and middle pitches in front of bottle cavities for the low ones.",
           (DEPTH, ("front", "float", (0.2, 0.8)), ("outer", "float", (0.05, 1.0)), ("inner", "float", (0.05, 1.0)),
            ("hole", "float", (2.0, 20.0)), ("plate", "float", (2.0, 10.0)), ("porosity", "float", (0.005, 0.2)),
            ("porous", "porous", None), ("rigid", "rigid", None)),
           hybrid),
)}


def decode(family, vector, brief):
    porous = brief.porous or list(FELTS)
    rigid = brief.rigid or list(RIGIDS)
    p = {}
    for u, (name, kind, spec) in zip(vector, family.genes):
        u = float(np.clip(u, 0, 1))
        if kind == "porous":
            p[name] = porous[min(int(u * len(porous)), len(porous) - 1)]
        elif kind == "rigid":
            p[name] = rigid[min(int(u * len(rigid)), len(rigid) - 1)]
        elif kind == "depth":
            lo, hi = spec
            p[name] = round((lo + u * (hi - lo)) * brief.max_depth, 1)
        else:
            lo, hi = spec
            p[name] = round(lo + u * (hi - lo), 4)
    return p


@dataclass
class FormResult:
    family: Family
    params: dict
    alpha: np.ndarray
    mass: float
    score: float
    fitness: float
    extra: dict = field(default_factory=dict)


def assess(family, params, brief, weights=None):
    weights = band_weights(brief.target) if weights is None else weights
    alpha, mass, extra = family.model(params, brief)
    score = float(np.sum(weights * alpha))
    fitness = score - (max(0.0, mass / brief.mass_budget - 1) + 0.1 * (mass > brief.mass_budget))
    return FormResult(family, params, alpha, float(mass), score, fitness, extra)


def breed(family, brief, population=30, generations=25, seed=1):
    weights = band_weights(brief.target)
    archive = {}

    def fitness(vector):
        params = decode(family, vector, brief)
        key = tuple(sorted(params.items()))
        if key not in archive:
            archive[key] = assess(family, params, brief, weights)
        return archive[key].fitness

    history = evolve(fitness, len(family.genes), population, generations, seed)
    return max(archive.values(), key=lambda r: r.fitness), history


def compare(brief, population=30, generations=25, seed=1, progress=None):
    """Best design of every family, ranked by fitness."""
    keys = brief.families or list(FAMILIES)
    results = []
    for key in keys:
        best, _ = breed(FAMILIES[key], brief, population, generations, seed)
        results.append(best)
        if progress:
            progress(best)
    return sorted(results, key=lambda r: r.fitness, reverse=True)


def material_matrix(brief, family_keys=("spikes", "sponge", "cloud", "hybrid"), population=20, generations=15):
    """Best score of each porous material inside each family: which material suits which shape."""
    table = {}
    for m in brief.porous or list(FELTS):
        sub = FormsBrief(brief.max_depth, brief.mass_budget, brief.target, [m], brief.rigid)
        table[m] = {k: breed(FAMILIES[k], sub, population, generations)[0] for k in family_keys}
    return table
