"""Console interface for running and verifying example thermal models."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from . import Boundary, Material, Slab, SolverOptions, simulate, isothermal_wall_reference
from .config import load_scenario
from .plotting import profile_plot, energy_plot, history_plot, verification_plot, convergence_plot


def run_case(config: Path, output: Path, no_plots: bool = False) -> dict:
    title, slab, options = load_scenario(config)
    output.mkdir(parents=True, exist_ok=True)
    result = simulate(slab, options)
    with (output / "temperatures.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "depth_m", "temperature_c"])
        for t, profile in zip(result.times_s, result.temperatures_c):
            for x, temp in zip(result.positions_m, profile):
                writer.writerow([f"{t:.10g}", f"{x:.10g}", f"{temp:.12g}"])

    np.savez_compressed(output / "simulation.npz", x_m=result.positions_m, times_s=result.times_s,
                        temperatures_c=result.temperatures_c, energy_j_m2=result.stored_energy_j_m2,
                        net_heat_j_m2=result.cumulative_net_heat_j_m2, balance_j_m2=result.energy_balance_error_j_m2)
    summary = {
        "title": title,
        "method": options.method,
        "material": slab.material.name,
        "diffusivity_m2_s": slab.material.diffusivity_m2_s,
        "thickness_m": slab.thickness_m,
        "cells": slab.cells,
        "dt_s": options.dt_s,
        "duration_s": float(result.times_s[-1]),
        "final_min_temperature_c": float(result.final_temperature_c.min()),
        "final_max_temperature_c": float(result.final_temperature_c.max()),
        "final_mean_temperature_c": float(result.final_temperature_c.mean()),
        "max_abs_energy_balance_error_j_m2": result.max_abs_balance_error_j_m2,
        "data_provenance": "Independent model; illustrative generic materials; no mission data",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if not no_plots:
        profile_plot(result, output / "temperature_profiles.png", title)
        history_plot(result, output / "temperature_history.png")
        energy_plot(result, output / "energy_balance.png")
    return summary


def verify_case(output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    material = Material(45.0, 7800.0, 500.0, "Illustrative generic alloy")
    cells_list = [20, 40, 80, 160]
    errors = []
    selected = None
    for cells in cells_list:
        slab = Slab(0.02, cells, material, Boundary.fixed(100.0), Boundary.fixed(100.0), 20.0)
        r = simulate(slab, SolverOptions(dt_s=0.025, steps=160, save_every=160))
        exact = isothermal_wall_reference(r.positions_m, r.times_s[-1], slab.thickness_m,
                                           material.diffusivity_m2_s, 20.0, 100.0)
        rmse = float(np.sqrt(np.mean((exact - r.final_temperature_c)**2)))
        errors.append(rmse)
        if cells == 80:
            selected = (r.positions_m, r.final_temperature_c, exact)
    verification_plot(*selected, output / "analytical_comparison.png")
    convergence_plot(cells_list, errors, output / "grid_convergence.png")
    report = {
        "reference": "Homogeneous 1D slab; equal fixed-temperature surfaces; Fourier-series solution",
        "duration_s": 4.0,
        "dt_s": 0.025,
        "cells": cells_list,
        "rmse_c": errors,
        "observed_convergence_ratio_20_to_40": errors[0] / errors[1],
        "observed_convergence_ratio_40_to_80": errors[1] / errors[2],
        "observed_convergence_ratio_80_to_160": errors[2] / errors[3],
    }
    (output / "verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="thermalphysics", description="Generic 1D transient heat-conduction workbench")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run a JSON scenario")
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--output", type=Path, default=Path("results/simulation"))
    run.add_argument("--no-plots", action="store_true")
    verify = sub.add_parser("verify", help="Compare against an analytical Fourier-series solution")
    verify.add_argument("--output", type=Path, default=Path("results/verification"))
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            report = run_case(args.config, args.output, args.no_plots)
        else:
            report = verify_case(args.output)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(report, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
