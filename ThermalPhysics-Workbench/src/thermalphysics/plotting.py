"""Headless scientific graphics for reproducible reports."""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from .solver import SimulationResult


def _save(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=160, bbox_inches="tight")
    plt.close()


def profile_plot(result: SimulationResult, path: str | Path, title: str = "Transient conduction") -> None:
    path = Path(path)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    take = np.unique(np.linspace(0, len(result.times_s)-1, min(6, len(result.times_s))).round().astype(int))
    for idx in take:
        ax.plot(result.positions_m * 1e3, result.temperatures_c[idx], label=f"t = {result.times_s[idx]:g} s", lw=2)
    ax.set(xlabel="Depth (mm)", ylabel="Temperature (°C)", title=title)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    _save(path)


def energy_plot(result: SimulationResult, path: str | Path) -> None:
    path = Path(path)
    delta = result.stored_energy_j_m2 - result.stored_energy_j_m2[0]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(result.times_s, delta, label="Δ internal energy", lw=2)
    ax.plot(result.times_s, result.cumulative_net_heat_j_m2, "--", label="Integrated net heat", lw=2)
    ax.set(xlabel="Time (s)", ylabel="Energy per unit area (J/m²)", title="Discrete energy-balance check")
    ax.legend()
    ax.grid(alpha=0.25)
    _save(path)


def history_plot(result: SimulationResult, path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for cell, label in [(0, "Left cell"), (result.positions_m.size//2, "Midplane"), (-1, "Right cell")]:
        ax.plot(result.times_s, result.temperatures_c[:, cell], lw=2, label=label)
    ax.set(xlabel="Time (s)", ylabel="Temperature (°C)", title="Selected cell temperature histories")
    ax.grid(alpha=0.25)
    ax.legend()
    _save(Path(path))


def verification_plot(x: np.ndarray, numerical: np.ndarray, analytical: np.ndarray, path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x * 1e3, analytical, lw=2, label="Analytical Fourier series")
    ax.plot(x * 1e3, numerical, "--", lw=2, label="Finite-volume Crank–Nicolson")
    ax.set(xlabel="Depth (mm)", ylabel="Temperature (°C)", title="Verification against an analytical solution")
    ax.grid(alpha=0.25)
    ax.legend()
    _save(Path(path))


def convergence_plot(cells: list[int], errors: list[float], path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.loglog(cells, errors, "o-", lw=2, label="RMSE at fixed time step")
    ax.loglog(cells, errors[0] * (cells[0] / np.asarray(cells))**2, "--", label="Reference: O(Δx²)")
    ax.set(xlabel="Number of cells", ylabel="RMSE (°C)", title="Spatial grid convergence")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    _save(Path(path))
