from pathlib import Path
import json
from thermalphysics.cli import run_case, verify_case

ROOT = Path(__file__).resolve().parents[1]


def test_run_produces_artifacts(tmp_path):
    summary = run_case(ROOT / "examples/insulated_internal_heating.json", tmp_path)
    assert summary["duration_s"] == 40
    assert (tmp_path / "temperatures.csv").exists()
    assert (tmp_path / "summary.json").exists()
    assert (tmp_path / "temperature_profiles.png").exists()
    assert (tmp_path / "energy_balance.png").exists()
    assert (tmp_path / "simulation.npz").exists()
    assert json.loads((tmp_path / "summary.json").read_text())["cells"] == 50


def test_validation_reduces_rmse_with_grid_refinement(tmp_path):
    report = verify_case(tmp_path)
    assert report["rmse_c"][0] > report["rmse_c"][1] > report["rmse_c"][2] > report["rmse_c"][3]
    assert (tmp_path / "grid_convergence.png").exists()
    assert (tmp_path / "analytical_comparison.png").exists()
