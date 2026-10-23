"""Layered media for general absorber forms.

Every form is reduced to a stack of layers along the depth axis, each one an
equivalent fluid (complex density and bulk modulus):
- porous felt (Miki model);
- felt mixed side by side with air, for tapered or open shapes
  (1/rho and 1/K mixed by volume fraction: spikes, sponges, lattices);
- narrow slits (exact viscous and thermal losses), for coiled channels;
- perforated plates (Maa), for resonator necks.
Same exp(+j*omega*t) convention as acoustics.py.
"""
from __future__ import annotations

import numpy as np

from .acoustics import C0, RHO0, Z0, miki

ETA = 1.81e-5  # dynamic viscosity of air, Pa*s
PRANDTL = 0.71
GAMMA = 1.4
K0 = RHO0 * C0 ** 2  # adiabatic bulk modulus of air


def porous_medium(f, sigma, frame_density=None):
    """Equivalent fluid of a felt. With `frame_density` (kg/m3) the fibres may move
    with the sound (limp model, Panneton 2007), as in light hanging objects."""
    omega = 2 * np.pi * f
    zc, kc = miki(f, sigma)
    rho, k = zc * kc / omega, zc * omega / kc
    if frame_density is not None:
        rho_t = frame_density + RHO0
        rho = (rho_t * rho - RHO0 ** 2) / (rho_t + rho - 2 * RHO0)
    return rho, k


def mixed_medium(f, sigma, solid_fraction, frame_density=None):
    """Felt occupying `solid_fraction` of the cross-section, open air the rest."""
    s = float(np.clip(solid_fraction, 0.0, 1.0))
    rho_p, k_p = porous_medium(f, sigma, frame_density)
    rho = 1 / (s / rho_p + (1 - s) / RHO0)
    k = 1 / (s / k_p + (1 - s) / K0)
    return rho, k


def slit_medium(f, width):
    """Air in a slit of `width` metres, with boundary-layer losses."""
    omega = 2 * np.pi * f
    h = width / 2
    a = h * np.sqrt(1j * omega * RHO0 / ETA)
    b = h * np.sqrt(1j * omega * RHO0 * PRANDTL / ETA)
    rho = RHO0 / (1 - np.tanh(a) / a)
    k = GAMMA * 101325.0 / (1 + (GAMMA - 1) * np.tanh(b) / b)
    return rho, k


def perforated_plate(f, hole_diameter, thickness, porosity):
    """Series impedance of a perforated plate (Maa), normal incidence."""
    omega = 2 * np.pi * f
    kp = hole_diameter / 2 * np.sqrt(omega * RHO0 / ETA)
    r = (32 * ETA * thickness / (porosity * RHO0 * C0 * hole_diameter ** 2)
         * (np.sqrt(1 + kp ** 2 / 32) + np.sqrt(2) / 32 * kp * hole_diameter / thickness))
    x = omega * thickness / (porosity * C0) * (1 + 1 / np.sqrt(9 + kp ** 2 / 2) + 0.85 * hole_diameter / thickness)
    return Z0 * (r + 1j * x)


def layer(f, theta, rho, k, d):
    """Transfer matrix (pressure, normal velocity) of one layer at oblique incidence."""
    omega = 2 * np.pi * f
    kx = omega / C0 * np.sin(theta)
    kz = np.sqrt(omega ** 2 * rho / k - kx ** 2 + 0j)
    kz = np.where(kz.imag > 0, -kz, kz)
    zz = omega * rho / kz
    c, s = np.cos(kz * d), np.sin(kz * d)
    return c, 1j * zz * s, 1j * s / zz, c


def series(z):
    """Transfer matrix of a thin series impedance (perforated plate)."""
    one = np.ones_like(z)
    return one, z, 0 * z, one


def multiply(a, b):
    return (a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3])


def stack(matrices):
    total = matrices[0]
    for m in matrices[1:]:
        total = multiply(total, m)
    return total


def admittance_on(matrix, y_back):
    """Surface admittance of a stack (front to back) closed by admittance `y_back`."""
    m11, m12, m21, m22 = matrix
    return (m21 + m22 * y_back) / (m11 + m12 * y_back)


def alpha_from_admittance(y, theta):
    cos_t = np.cos(theta)
    r = (cos_t - Z0 * y) / (cos_t + Z0 * y)
    return np.clip(1 - np.abs(r) ** 2, 0, 1)


def alpha_free(matrix, theta):
    """Absorbed share of sound hitting a layer with open air on both sides (1 - R^2 - T^2)."""
    m11, m12, m21, m22 = matrix
    z0z = Z0 / np.cos(theta)
    denom = m11 + m12 / z0z + z0z * m21 + m22
    r = (m11 + m12 / z0z - z0z * m21 - m22) / denom
    t = 2 / denom
    return np.clip(1 - np.abs(r) ** 2 - np.abs(t) ** 2, 0, 1)


def paris(alpha, theta):
    """Random-incidence average over angles (rows = frequencies, columns = angles)."""
    w = np.sin(2 * theta)
    return np.trapezoid(alpha * w, theta, axis=-1) / np.trapezoid(w, theta)
