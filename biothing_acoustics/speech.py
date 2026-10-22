"""Speech privacy across a desk screen.

Direct and diffracted paths only (over the top, around both run ends and
through the panel). Reflections from ceiling, walls and the desk surface are
not modelled, so in low rooms with a hard ceiling real performance is lower.
"""
from __future__ import annotations

import math

import numpy as np

from .acoustics import C0

OCTAVES = np.array([125, 250, 500, 1000, 2000, 4000, 8000], dtype=float)
# Normal-effort speech, free field at 1 m (ISO 3382-3, table 1).
SPEECH_LEVEL_1M = np.array([49.9, 54.3, 58.0, 52.0, 44.8, 38.8, 33.5])
# Octave-band importance for the speech intelligibility index (ANSI S3.5).
BAND_IMPORTANCE = np.array([0.0, 0.0617, 0.1671, 0.2373, 0.2648, 0.2142, 0.0549])
A_WEIGHTING = np.array([-16.1, -8.6, -3.2, 0.0, 1.2, 1.0, -1.1])
MAX_BARRIER_ATTENUATION = 20.0  # dB, practical limit for a thin screen


def barrier_attenuation(n):
    """Single-edge diffraction attenuation from the Fresnel number (Kurze-Anderson)."""
    n = np.asarray(n, dtype=float)
    out = np.zeros_like(n)
    pos = n > 0
    x = np.sqrt(2 * np.pi * n[pos])
    out[pos] = 5 + 20 * np.log10(x / np.tanh(x))
    out[n == 0] = 5.0
    shadow_edge = (n < 0) & (n > -0.2)  # line of sight just clears the edge
    x = np.sqrt(2 * np.pi * -n[shadow_edge])
    out[shadow_edge] = 5 + 20 * np.log10(x / np.tan(x))
    return np.clip(out, 0, MAX_BARRIER_ATTENUATION)


def top_path_difference(a, b, hs, hr, top):
    """Signed path difference over the screen top (negative = line of sight open)."""
    direct = math.hypot(a + b, hr - hs)
    over = math.hypot(a, top - hs) + math.hypot(b, top - hr)
    sight = hs + (hr - hs) * a / (a + b)
    return (over - direct) * (1 if top >= sight else -1)


def screen_attenuation(scenario, top_mm, run_length_mm, tl, freqs=OCTAVES):
    """Attenuation (dB per band) of the direct speech path caused by a screen run."""
    a = scenario.talker_distance / 1000
    b = scenario.listener_distance / 1000
    hs = scenario.mouth_height / 1000
    hr = scenario.ear_height / 1000
    wavelength = C0 / np.asarray(freqs, dtype=float)
    n_top = 2 * top_path_difference(a, b, hs, hr, top_mm / 1000) / wavelength
    half = run_length_mm / 2000
    side = math.hypot(a, half) + math.hypot(b, half) - (a + b)
    n_side = 2 * side / wavelength
    tau = (10 ** (-barrier_attenuation(n_top) / 10)
           + 2 * 10 ** (-barrier_attenuation(n_side) / 10)
           + 10 ** (-np.asarray(tl) / 10))
    return -10 * np.log10(tau)


def sii(levels, noise):
    """Simplified speech intelligibility index from band signal-to-noise ratios."""
    audible = np.clip((np.asarray(levels) - np.asarray(noise) + 15) / 30, 0, 1)
    return float(np.sum(BAND_IMPORTANCE * audible))


def a_weighted(levels):
    return float(10 * np.log10(np.sum(10 ** ((np.asarray(levels) + A_WEIGHTING) / 10))))


def privacy(scenario, attenuation, noise):
    """Speech level and intelligibility at the listener, without and with the screen."""
    r = math.hypot((scenario.talker_distance + scenario.listener_distance) / 1000,
                   (scenario.ear_height - scenario.mouth_height) / 1000)
    open_levels = SPEECH_LEVEL_1M - 20 * math.log10(r)
    screened = open_levels - attenuation
    sii_open, sii_screen = sii(open_levels, noise), sii(screened, noise)
    return {
        "speech_open_dba": a_weighted(open_levels),
        "speech_screened_dba": a_weighted(screened),
        "speech_reduction_dba": a_weighted(open_levels) - a_weighted(screened),
        "sii_open": sii_open,
        "sii_screened": sii_screen,
        "privacy_gain": (sii_open - sii_screen) / sii_open if sii_open > 0 else 0.0,
        "levels_open": open_levels,
        "levels_screened": screened,
    }
