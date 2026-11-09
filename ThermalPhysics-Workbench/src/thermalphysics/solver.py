"""Conservative cell-centered finite-volume heat solver.

All lengths in m, time in s, temperatures in degC, fluxes in W/m^2.
No mission-specific properties, geometry, or data are included.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import solve_banded

from .model import Boundary, Slab, SolverOptions


@dataclass
class SimulationResult:
    positions_m: np.ndarray
    times_s: np.ndarray
    temperatures_c: np.ndarray  # (snapshots, cells)
    stored_energy_j_m2: np.ndarray
    cumulative_net_heat_j_m2: np.ndarray
    energy_balance_error_j_m2: np.ndarray
    method: str

    @property
    def final_temperature_c(self) -> np.ndarray:
        return self.temperatures_c[-1]

    @property
    def max_abs_balance_error_j_m2(self) -> float:
        return float(np.max(np.abs(self.energy_balance_error_j_m2)))


def _surface_conductance(boundary: Boundary, conductivity: float, dx: float) -> float:
    """Conductance from cell center to prescribed boundary temperature, W/m^2/K."""
    if boundary.kind == "fixed":
        return 2.0 * conductivity / dx
    if boundary.kind == "convection":
        return 1.0 / (dx / (2.0 * conductivity) + 1.0 / boundary.h_w_m2k)
    return 0.0


def _operator(slab: Slab) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute A and b of dT/dt=A@T+b, storing only tridiagonal entries."""
    n, dx = slab.cells, slab.dx_m
    k = slab.material.conductivity_w_mk
    cap = slab.material.density_kg_m3 * slab.material.heat_capacity_j_kgk * dx
    face_conductance = k / dx
    rate = face_conductance / cap

    diag = np.full(n, -2.0 * rate)
    diag[0] = -rate  # one neighboring cell + exterior boundary
    diag[-1] = -rate
    lower = np.full(n - 1, rate)
    upper = np.full(n - 1, rate)
    source = np.full(n, slab.heat_source_w_m3 * dx / cap)

    for index, boundary in ((0, slab.left), (-1, slab.right)):
        g_surface = _surface_conductance(boundary, k, dx)
        if boundary.kind in ("fixed", "convection"):
            ambient = boundary.temperature_c if boundary.kind == "fixed" else boundary.ambient_c
            diag[index] -= g_surface / cap
            source[index] += g_surface * ambient / cap
        elif boundary.kind == "flux":
            source[index] -= boundary.outward_flux_w_m2 / cap
    return lower, diag, upper, source


def _boundary_outward_flux(boundary: Boundary, cell_temperature: float, k: float, dx: float) -> float:
    """Instantaneous outward heat flux at a slab boundary (W/m^2)."""
    if boundary.kind == "insulated":
        return 0.0
    if boundary.kind == "flux":
        return float(boundary.outward_flux_w_m2)
    g = _surface_conductance(boundary, k, dx)
    ref = boundary.temperature_c if boundary.kind == "fixed" else boundary.ambient_c
    return g * (cell_temperature - ref)


def net_heat_rate_w_m2(slab: Slab, temperature_c: np.ndarray) -> float:
    """Volumetric generation minus heat leaving the left and right surfaces."""
    k = slab.material.conductivity_w_mk
    left = _boundary_outward_flux(slab.left, temperature_c[0], k, slab.dx_m)
    right = _boundary_outward_flux(slab.right, temperature_c[-1], k, slab.dx_m)
    return slab.heat_source_w_m3 * slab.thickness_m - left - right


def _banded(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray) -> np.ndarray:
    n = diag.size
    ab = np.zeros((3, n), dtype=float)
    ab[0, 1:] = upper
    ab[1, :] = diag
    ab[2, :-1] = lower
    return ab


def _apply_operator(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = diag * x
    y[1:] += lower * x[:-1]
    y[:-1] += upper * x[1:]
    return y


def simulate(slab: Slab, options: SolverOptions) -> SimulationResult:
    """Run a conservative transient 1D heat-conduction simulation.

    Crank-Nicolson is A-stable for the linear model. Explicit Euler is
    restricted by the nonnegative-update sufficient time-step bound.
    """
    lower, diag, upper, source = _operator(slab)
    dt = options.dt_s
    if options.method == "explicit-euler":
        dt_max = 1.0 / np.max(-diag)
        if dt > dt_max * (1.0 + 1e-12):
            raise ValueError(
                f"Unstable explicit step dt={dt:g}s exceeds conservative "
                f"positivity limit {dt_max:.6g}s; lower dt or use crank-nicolson"
            )
    else:
        left_mat = _banded(-0.5 * dt * lower, 1.0 - 0.5 * dt * diag, -0.5 * dt * upper)

    temperature = np.full(slab.cells, slab.initial_temperature_c, dtype=float)
    positions = (np.arange(slab.cells) + 0.5) * slab.dx_m
    area_capacity = slab.material.density_kg_m3 * slab.material.heat_capacity_j_kgk * slab.dx_m
    start_energy = area_capacity * float(temperature.sum())
    cumulative = 0.0
    initial_net = net_heat_rate_w_m2(slab, temperature)

    times = [0.0]
    temperatures = [temperature.copy()]
    energies = [start_energy]
    cumulative_heats = [0.0]
    balances = [0.0]

    for step in range(1, options.steps + 1):
        if options.method == "crank-nicolson":
            rhs = temperature + 0.5 * dt * _apply_operator(lower, diag, upper, temperature) + dt * source
            newer = solve_banded((1, 1), left_mat, rhs, check_finite=False)
        else:
            newer = temperature + dt * (_apply_operator(lower, diag, upper, temperature) + source)

        if not np.all(np.isfinite(newer)):
            raise FloatingPointError("Simulation diverged to non-finite temperatures")
        new_net = net_heat_rate_w_m2(slab, newer)
        if options.method == "crank-nicolson":
            cumulative += 0.5 * dt * (initial_net + new_net)
        else:
            cumulative += dt * initial_net
        temperature = newer
        initial_net = new_net

        if step % options.save_every == 0 or step == options.steps:
            energy = area_capacity * float(temperature.sum())
            times.append(step * dt)
            temperatures.append(temperature.copy())
            energies.append(energy)
            cumulative_heats.append(cumulative)
            balances.append((energy - start_energy) - cumulative)

    return SimulationResult(
        positions_m=positions,
        times_s=np.asarray(times),
        temperatures_c=np.stack(temperatures),
        stored_energy_j_m2=np.asarray(energies),
        cumulative_net_heat_j_m2=np.asarray(cumulative_heats),
        energy_balance_error_j_m2=np.asarray(balances),
        method=options.method,
    )


def steady_state(slab: Slab) -> np.ndarray:
    """Solve A@T + b = 0. Requires at least one temperature-anchoring boundary."""
    if slab.left.kind not in {"fixed", "convection"} and slab.right.kind not in {"fixed", "convection"}:
        raise ValueError("Steady state is not uniquely defined without a fixed/convection boundary")
    lo, di, up, src = _operator(slab)
    return solve_banded((1, 1), _banded(lo, di, up), -src)
