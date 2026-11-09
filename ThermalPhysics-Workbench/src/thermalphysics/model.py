"""Validated, unit-explicit models for a generic one-dimensional slab."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite


def _finite(name: str, value: float) -> None:
    if not isinstance(value, (int, float)) or not isfinite(value):
        raise ValueError(f"{name} must be a finite number")


@dataclass(frozen=True)
class Material:
    """Homogeneous material with temperature-independent thermal properties (SI)."""

    conductivity_w_mk: float
    density_kg_m3: float
    heat_capacity_j_kgk: float
    name: str = "Generic material"

    def __post_init__(self) -> None:
        for key in ("conductivity_w_mk", "density_kg_m3", "heat_capacity_j_kgk"):
            value = getattr(self, key)
            _finite(key, value)
            if value <= 0:
                raise ValueError(f"{key} must be greater than zero")

    @property
    def diffusivity_m2_s(self) -> float:
        return self.conductivity_w_mk / (self.density_kg_m3 * self.heat_capacity_j_kgk)


@dataclass(frozen=True)
class Boundary:
    """Generic surface BC. Outward heat flux is positive (W/m^2)."""

    kind: str
    temperature_c: float | None = None
    h_w_m2k: float | None = None
    ambient_c: float | None = None
    outward_flux_w_m2: float | None = None

    def __post_init__(self) -> None:
        if self.kind not in {"fixed", "insulated", "convection", "flux"}:
            raise ValueError("Boundary kind must be fixed, insulated, convection, or flux")
        required = {
            "fixed": {"temperature_c"},
            "insulated": set(),
            "convection": {"h_w_m2k", "ambient_c"},
            "flux": {"outward_flux_w_m2"},
        }[self.kind]
        fields = {"temperature_c", "h_w_m2k", "ambient_c", "outward_flux_w_m2"}
        for key in fields:
            value = getattr(self, key)
            if key in required and value is None:
                raise ValueError(f"{self.kind} boundary requires {key}")
            if key not in required and value is not None:
                raise ValueError(f"{key} does not apply to {self.kind} boundary")
            if value is not None:
                _finite(key, value)
        if self.kind == "convection" and self.h_w_m2k <= 0:
            raise ValueError("Convection h_w_m2k must be positive")

    @classmethod
    def fixed(cls, temperature_c: float) -> "Boundary":
        return cls("fixed", temperature_c=temperature_c)

    @classmethod
    def insulated(cls) -> "Boundary":
        return cls("insulated")

    @classmethod
    def convection(cls, h_w_m2k: float, ambient_c: float) -> "Boundary":
        return cls("convection", h_w_m2k=h_w_m2k, ambient_c=ambient_c)

    @classmethod
    def flux(cls, outward_flux_w_m2: float) -> "Boundary":
        return cls("flux", outward_flux_w_m2=outward_flux_w_m2)


@dataclass(frozen=True)
class Slab:
    thickness_m: float
    cells: int
    material: Material
    left: Boundary
    right: Boundary
    initial_temperature_c: float
    heat_source_w_m3: float = 0.0

    def __post_init__(self) -> None:
        _finite("thickness_m", self.thickness_m)
        _finite("initial_temperature_c", self.initial_temperature_c)
        _finite("heat_source_w_m3", self.heat_source_w_m3)
        if self.thickness_m <= 0:
            raise ValueError("thickness_m must be positive")
        if isinstance(self.cells, bool) or not isinstance(self.cells, int) or self.cells < 3:
            raise ValueError("cells must be an integer >= 3")

    @property
    def dx_m(self) -> float:
        return self.thickness_m / self.cells


@dataclass(frozen=True)
class SolverOptions:
    dt_s: float
    steps: int
    save_every: int = 1
    method: str = "crank-nicolson"

    def __post_init__(self) -> None:
        _finite("dt_s", self.dt_s)
        if self.dt_s <= 0:
            raise ValueError("dt_s must be positive")
        if isinstance(self.steps, bool) or not isinstance(self.steps, int) or self.steps < 1:
            raise ValueError("steps must be a positive integer")
        if isinstance(self.save_every, bool) or not isinstance(self.save_every, int) or self.save_every < 1:
            raise ValueError("save_every must be a positive integer")
        if self.method not in {"crank-nicolson", "explicit-euler"}:
            raise ValueError("method must be crank-nicolson or explicit-euler")
