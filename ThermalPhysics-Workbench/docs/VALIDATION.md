# Verification and test strategy

**Verification** asks whether the equations are solved correctly. **Validation** asks whether the equations represent the real-world physics sufficiently well. This repository includes numerical **verification**, but does not claim experimental validation.

## Automated checks

`python -m pytest` covers:

1. Constructor input validation and nonphysical parameter rejection.
2. Preservation of a uniform initial field between insulated walls.
3. Exact energy gain for a slab with constant internal heating and insulated surfaces.
4. Cooling under two convective boundary conditions with a discrete energy-balance check.
5. Sign convention of an imposed outward surface heat flux.
6. Analytical Fourier-series comparison for a uniformly initialized slab with equal fixed-temperature walls.
7. Steady linear temperature profile between different prescribed wall temperatures.
8. Explicit Euler time-step rejection when the positivity bound is exceeded.
9. Convergence of explicit and implicit methods for short stable time steps.
10. CLI export of structured CSV, NPZ, JSON and plots.
11. Grid refinement reduces RMS error against the analytical reference.

## Reproduce analytical checks

```bash
thermalphysics verify --output results/verification
```

Generates:
- `verification.json` — RMS errors and convergence ratios
- `analytical_comparison.png` — numerical solution vs Fourier series
- `grid_convergence.png` — spatial error vs cell count

Input data are generic **illustrative material parameters**, not empirical certification-grade data.

## Reproduced numerical results (Python 3.13, release 1.0.0)

The included run reports the following root-mean-square errors against the analytical Fourier-series solution:

| Control volumes | RMSE (°C) | Error reduction from previous grid |
|---:|---:|---:|
| 20 | 0.0776943 | — |
| 40 | 0.0193233 | 4.02× |
| 80 | 0.00474658 | 4.07× |
| 160 | 0.00117371 | 4.04× |

Approximately fourfold reductions for each doubling of cell count demonstrate the expected second-order **spatial** convergence in this benchmark. Reproduce rather than assume these results on new environments. Raw metrics are in [`verification-report.json`](verification-report.json).

## Interpretation

For fixed $Δt$, second-order spatial convergence should be seen while spatial error dominates. Eventually temporal error and floating-point limits prevent further improvement. Absolute errors depend on time step, wall temperatures and chosen physical parameters.

Crank–Nicolson is second order in time for this linear system but large time steps still lose precision. No plot should be interpreted as proof that a heat-protection system is flightworthy.
