import numpy as np
import pytest
from thermalphysics import (
    Material, Boundary, Slab, SolverOptions, simulate, steady_state,
    isothermal_wall_reference, net_heat_rate_w_m2
)

M = Material(45, 7800, 500, "Illustrative solid")


def slab(**changes):
    d = dict(thickness_m=.02, cells=60, material=M, left=Boundary.insulated(),
             right=Boundary.insulated(), initial_temperature_c=20.0, heat_source_w_m3=0)
    d.update(changes)
    return Slab(**d)


def test_insulated_uniform_stays_constant():
    r = simulate(slab(), SolverOptions(.2, 100, 10))
    np.testing.assert_allclose(r.temperatures_c, 20, rtol=0, atol=1e-10)
    assert r.max_abs_balance_error_j_m2 < 1e-5


def test_internal_generation_energy_and_uniform_warming():
    s = slab(heat_source_w_m3=125000)
    r = simulate(s, SolverOptions(.2, 100, 5))
    expected = 20 + 125000 * 20 / (M.density_kg_m3 * M.heat_capacity_j_kgk)
    np.testing.assert_allclose(r.final_temperature_c, expected, atol=1e-8)
    assert r.max_abs_balance_error_j_m2 < 1e-5


def test_fixed_wall_matches_analytical_solution():
    s = slab(cells=100, left=Boundary.fixed(100), right=Boundary.fixed(100))
    r = simulate(s, SolverOptions(.025, 160, 160))
    exact = isothermal_wall_reference(r.positions_m, 4.0, .02, M.diffusivity_m2_s, 20, 100)
    assert np.sqrt(np.mean((r.final_temperature_c-exact)**2)) < 0.02
    assert r.max_abs_balance_error_j_m2 < 1e-4


def test_convection_cooling_and_energy_balance():
    s = slab(initial_temperature_c=180, left=Boundary.convection(80, 20),
             right=Boundary.convection(125, 20))
    r = simulate(s, SolverOptions(.4, 200, 10))
    assert r.final_temperature_c.mean() < 180
    assert r.final_temperature_c.min() > 19.99
    assert r.max_abs_balance_error_j_m2 < 1e-4


def test_outward_flux_removes_energy():
    s = slab(left=Boundary.flux(1000))
    assert net_heat_rate_w_m2(s, np.full(s.cells,20)) == -1000
    r = simulate(s, SolverOptions(.1, 100, 10))
    expected_drop = 1000 * 10 / (M.density_kg_m3 * M.heat_capacity_j_kgk * s.thickness_m)
    assert r.final_temperature_c.mean() == pytest.approx(20-expected_drop, abs=1e-8)
    assert r.max_abs_balance_error_j_m2 < 1e-5


def test_explicit_stability_guard():
    with pytest.raises(ValueError, match="Unstable explicit step"):
        simulate(slab(left=Boundary.fixed(60)), SolverOptions(100, 10, method="explicit-euler"))


def test_explicit_agrees_with_implicit_for_small_step():
    s = slab(cells=20, left=Boundary.fixed(90), right=Boundary.fixed(40))
    explicit = simulate(s, SolverOptions(.00005, 2000, 2000, "explicit-euler"))
    implicit = simulate(s, SolverOptions(.00005, 2000, 2000, "crank-nicolson"))
    np.testing.assert_allclose(explicit.final_temperature_c, implicit.final_temperature_c, atol=1e-2)


def test_steady_state_is_linear_with_two_fixed_walls():
    s = slab(left=Boundary.fixed(100), right=Boundary.fixed(40))
    t = steady_state(s)
    np.testing.assert_allclose(t, 100 + (40-100) * ((np.arange(s.cells)+.5)/s.cells), atol=1e-9)


def test_steady_state_requires_anchor():
    with pytest.raises(ValueError):
        steady_state(slab())


def test_analytic_initial_temperature():
    x = np.linspace(0.001, 0.019, 5)
    np.testing.assert_allclose(isothermal_wall_reference(x, 0, .02, 1e-5, 20, 100), 20)
