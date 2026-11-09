"""Textbook Fourier-series reference for a simple slab boundary-value problem."""
from __future__ import annotations
import numpy as np


def isothermal_wall_reference(
    positions_m: np.ndarray,
    time_s: float,
    thickness_m: float,
    diffusivity_m2_s: float,
    initial_temperature_c: float,
    wall_temperature_c: float,
    terms: int = 400,
) -> np.ndarray:
    """Exact series: initially uniform slab; both walls fixed at same T.

    Assumes constant thermal diffusivity and no source. At exactly t=0
    returns the interior initial condition; the boundary is discontinuous.
    """
    if time_s < 0 or thickness_m <= 0 or diffusivity_m2_s <= 0 or terms < 1:
        raise ValueError("Expect t>=0, thickness>0, diffusivity>0, terms>=1")
    x = np.asarray(positions_m, dtype=float)
    if np.any((x < 0) | (x > thickness_m)):
        raise ValueError("positions must be within the slab")
    if time_s == 0:
        return np.full_like(x, initial_temperature_c)
    odd = (2 * np.arange(terms) + 1).astype(float)
    p = odd[:, None] * np.pi / thickness_m
    coefficient = (4.0 / (odd * np.pi))[:, None]
    factor = np.exp(-diffusivity_m2_s * (odd * np.pi / thickness_m) ** 2 * time_s)[:, None]
    theta = np.sum(coefficient * np.sin(p * x[None, :]) * factor, axis=0)
    return wall_temperature_c + (initial_temperature_c - wall_temperature_c) * theta
