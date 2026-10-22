"""Sound absorption and sound reduction models.

Convention: time dependence exp(+j*omega*t), so lossy wavenumbers have a
negative imaginary part.

Absorption of a sculpted face is predicted for an infinite surface. Small
panels measured in a reverberation room (ISO 354) usually score higher
because of edge diffraction, so treat the numbers as conservative and
comparative, then verify with measurements.
"""
from __future__ import annotations

import numpy as np

RHO0 = 1.204  # kg/m3, air at 20 C
C0 = 343.2  # m/s
Z0 = RHO0 * C0

THIRD_OCTAVES = np.array([100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000,
                          1250, 1600, 2000, 2500, 3150, 4000, 5000], dtype=float)
RATING_OCTAVES = (250, 500, 1000, 2000, 4000)  # ISO 11654
COMPRESSION_EXPONENT = 1.53  # sigma ~ density^1.53 for fibrous materials (Bies & Hansen)


def miki(f, sigma):
    """Characteristic impedance and complex wavenumber of a porous layer (Miki, 1990)."""
    f = np.asarray(f, dtype=float)
    x = f / np.asarray(sigma, dtype=float)
    zc = Z0 * (1 + 0.070 * x ** -0.632 - 1j * 0.107 * x ** -0.632)
    kc = (2 * np.pi * f / C0) * (1 + 0.109 * x ** -0.618 - 1j * 0.160 * x ** -0.618)
    return zc, kc


def compressed_sigma(sigma0, ratio):
    """Flow resistivity of a felt pressed to `ratio` (0..1] of its original thickness."""
    return np.asarray(sigma0) * np.asarray(ratio, dtype=float) ** -COMPRESSION_EXPONENT


def surface_impedance(f, theta, thickness, gap, sigma):
    """Normal surface impedance of a felt layer over an air gap on a rigid core.

    All arguments broadcast together (SI units: Hz, rad, m, m, Pa*s/m2).
    The felt is modelled as extended-reacting, the gap as locally reacting.
    """
    k0 = 2 * np.pi * np.asarray(f, dtype=float) / C0
    cos_t = np.cos(theta)
    kx = k0 * np.sin(theta)
    # Admittance of the air gap closed by the rigid core (zero when gap == 0).
    yb = 1j * cos_t * np.tan(k0 * cos_t * gap) / Z0
    zc, kc = miki(f, sigma)
    kz = np.sqrt(kc ** 2 - kx ** 2)
    zcz = zc * kc / kz
    t = np.tan(kz * thickness)
    return zcz * (1 + 1j * zcz * t * yb) / (zcz * yb + 1j * t)


def panel_absorption(area_fraction, thickness, gap, sigma, freqs=THIRD_OCTAVES, n_theta=46):
    """Absorption of a face made of cells with different felt thickness, gap and sigma.

    Cells are combined by area-weighted admittance (valid while cells are small
    compared with the wavelength, i.e. up to a few kHz for 20-80 mm cells).
    Returns normal-incidence and random-incidence (Paris) absorption per frequency.
    """
    w = np.asarray(area_fraction, dtype=float)[:, None, None]
    d = np.asarray(thickness, dtype=float)[:, None, None]
    g = np.asarray(gap, dtype=float)[:, None, None]
    s = np.asarray(sigma, dtype=float)[:, None, None]
    f = np.asarray(freqs, dtype=float)[None, :, None]
    theta = np.linspace(0, np.radians(89), n_theta)
    zs = surface_impedance(f, theta[None, None, :], d, g, s)
    y = np.sum(w / zs, axis=0)  # (F, T)
    zeff = 1 / y
    cos_t = np.cos(theta)[None, :]
    r = (zeff * cos_t - Z0) / (zeff * cos_t + Z0)
    alpha = np.clip(1 - np.abs(r) ** 2, 0, 1)
    weight = np.sin(2 * theta)
    random = np.trapezoid(alpha * weight, theta, axis=1) / np.trapezoid(weight, theta)
    return {"freqs": np.asarray(freqs, dtype=float), "normal": alpha[:, 0], "random": random}


def practical_alpha(alpha_thirds, freqs=THIRD_OCTAVES):
    """Practical absorption coefficient per octave (ISO 11654): mean of three
    one-third octaves, rounded to 0.05 and capped at 1.00."""
    freqs = np.asarray(freqs)
    out = {}
    for octave in (125,) + RATING_OCTAVES:
        i = int(np.argmin(np.abs(freqs - octave)))
        mean = float(np.mean(alpha_thirds[i - 1:i + 2]))
        out[octave] = min(1.0, round(round(mean / 0.05) * 0.05, 2))
    return out


def weighted_alpha(alpha_p):
    """Weighted sound absorption coefficient, shape indicators and class (ISO 11654).

    `alpha_p` holds the practical coefficients at 250, 500, 1000, 2000 and 4000 Hz.
    """
    reference = (80, 100, 100, 100, 90)
    ap = [int(round(a * 100)) for a in alpha_p]
    for aw in range(100, -1, -5):
        shifted = [r - (100 - aw) for r in reference]
        if sum(max(0, r - a) for r, a in zip(shifted, ap)) <= 10:
            break
    shape = ""
    for letter, r, a in zip("LMMHH", shifted, ap):
        if a - r >= 25 and letter not in shape:
            shape += letter
    if aw >= 90:
        cls = "A"
    elif aw >= 80:
        cls = "B"
    elif aw >= 60:
        cls = "C"
    elif aw >= 30:
        cls = "D"
    elif aw >= 15:
        cls = "E"
    else:
        cls = "not classified"
    return aw / 100, shape, cls


def nrc(alpha_thirds, freqs=THIRD_OCTAVES):
    """Noise reduction coefficient: mean of 250-2000 Hz octaves, rounded to 0.05."""
    freqs = np.asarray(freqs)
    bands = []
    for octave in (250, 500, 1000, 2000):
        i = int(np.argmin(np.abs(freqs - octave)))
        bands.append(np.mean(alpha_thirds[i - 1:i + 2]))
    return round(round(float(np.mean(bands)) / 0.05) * 0.05, 2)


def mass_law_tl(surface_mass, f):
    """Field-incidence mass law transmission loss in dB (ignores coincidence)."""
    return np.maximum(0.0, 20 * np.log10(surface_mass * np.asarray(f, dtype=float)) - 47)


def composite_tl(tl, leak_fraction):
    """Transmission loss of a panel with open slits (joints, gap under the screen)."""
    tau = (1 - leak_fraction) * 10 ** (-np.asarray(tl) / 10) + leak_fraction
    return -10 * np.log10(tau)
