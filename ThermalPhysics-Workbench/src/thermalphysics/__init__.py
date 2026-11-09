"""ThermalPhysics Workbench: original, independent 1D heat-transfer teaching toolkit."""
from .model import Boundary, Material, Slab, SolverOptions
from .solver import SimulationResult, net_heat_rate_w_m2, simulate, steady_state
from .analytic import isothermal_wall_reference

__version__ = "1.0.0"
__all__ = ["Boundary", "Material", "Slab", "SolverOptions", "SimulationResult", "simulate", "steady_state", "net_heat_rate_w_m2", "isothermal_wall_reference"]
