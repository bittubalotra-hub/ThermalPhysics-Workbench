"""JSON scenario input with explicit, readable SI units."""
from __future__ import annotations
import json
from pathlib import Path
from .model import Material, Boundary, Slab, SolverOptions


def _read_boundary(data: dict) -> Boundary:
    kind = data["type"]
    if kind == "fixed":
        return Boundary.fixed(float(data["temperature_c"]))
    if kind == "insulated":
        return Boundary.insulated()
    if kind == "convection":
        return Boundary.convection(float(data["h_w_m2k"]), float(data["ambient_c"]))
    if kind == "flux":
        return Boundary.flux(float(data["outward_flux_w_m2"]))
    raise ValueError(f"Unknown boundary type: {kind}")


def load_scenario(path: str | Path) -> tuple[str, Slab, SolverOptions]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    m = raw["material"]
    g = raw["geometry"]
    b = raw["boundaries"]
    s = raw["solver"]
    material = Material(
        name=str(m.get("name", "Generic material")),
        conductivity_w_mk=float(m["conductivity_w_mk"]),
        density_kg_m3=float(m["density_kg_m3"]),
        heat_capacity_j_kgk=float(m["heat_capacity_j_kgk"]),
    )
    slab = Slab(
        thickness_m=float(g["thickness_m"]),
        cells=g["cells"],
        material=material,
        left=_read_boundary(b["left"]),
        right=_read_boundary(b["right"]),
        initial_temperature_c=float(raw["initial_temperature_c"]),
        heat_source_w_m3=float(raw.get("heat_source_w_m3", 0.0)),
    )
    options = SolverOptions(
        dt_s=float(s["dt_s"]),
        steps=s["steps"],
        save_every=s.get("save_every", 1),
        method=s.get("method", "crank-nicolson"),
    )
    return str(raw.get("title", Path(path).stem)), slab, options
