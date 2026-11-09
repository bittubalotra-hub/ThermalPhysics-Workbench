"""Independent cross-language check of C++ and Python finite-volume update."""
from pathlib import Path
import argparse
import csv
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from thermalphysics import Material, Boundary, Slab, SolverOptions, simulate  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description="Compare C++ CSV profile with Python solver")
    p.add_argument("csv", type=Path)
    p.add_argument("--cells", type=int, default=80)
    p.add_argument("--steps", type=int, default=8000)
    p.add_argument("--dt", type=float, default=0.0005)
    args = p.parse_args()
    with args.csv.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    x_cpp = np.array([float(row["depth_m"]) for row in rows])
    t_cpp = np.array([float(row["temperature_c"]) for row in rows])
    material = Material(45, 7800, 500, "Generic example material")
    slab = Slab(.02, args.cells, material, Boundary.fixed(100), Boundary.fixed(100), 20)
    result = simulate(slab, SolverOptions(args.dt, args.steps, args.steps, "explicit-euler"))
    np.testing.assert_allclose(x_cpp, result.positions_m, atol=1e-13, rtol=0)
    np.testing.assert_allclose(t_cpp, result.final_temperature_c, atol=1e-9, rtol=0)
    print(f"PASS: cross-language max temperature difference = {np.max(np.abs(t_cpp - result.final_temperature_c)):.3e} °C")

if __name__ == "__main__":
    main()
