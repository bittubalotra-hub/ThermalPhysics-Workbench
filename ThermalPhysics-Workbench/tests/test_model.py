import pytest
from thermalphysics import Material, Boundary, Slab, SolverOptions


def test_properties_and_diffusivity():
    m = Material(40, 8000, 500)
    assert m.diffusivity_m2_s == pytest.approx(1e-5)


@pytest.mark.parametrize("kwargs", [
    {"conductivity_w_mk": 0, "density_kg_m3": 8000, "heat_capacity_j_kgk": 500},
    {"conductivity_w_mk": 40, "density_kg_m3": -2, "heat_capacity_j_kgk": 500},
    {"conductivity_w_mk": 40, "density_kg_m3": 8000, "heat_capacity_j_kgk": float("nan")},
])
def test_bad_material(kwargs):
    with pytest.raises(ValueError):
        Material(**kwargs)


def test_boundary_validation():
    assert Boundary.insulated().kind == "insulated"
    with pytest.raises(ValueError):
        Boundary.convection(0, 20)
    with pytest.raises(ValueError):
        Boundary("fixed", temperature_c=20, ambient_c=20)


def test_options_validation():
    with pytest.raises(ValueError):
        SolverOptions(dt_s=0, steps=1)
    with pytest.raises(ValueError):
        SolverOptions(dt_s=1, steps=0)
    with pytest.raises(ValueError):
        SolverOptions(dt_s=1, steps=1, method="unknown")


def test_slab_validation():
    m = Material(40, 8000, 500)
    with pytest.raises(ValueError):
        Slab(0.01, 2, m, Boundary.insulated(), Boundary.insulated(), 20)
