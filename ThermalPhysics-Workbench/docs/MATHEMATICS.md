# Mathematical model and numerical formulation

The code solves the **one-dimensional, transient, linear heat equation** in a homogeneous slab of thickness $L$:

$$\rho c_p \frac{\partial T}{\partial t}=k\frac{\partial^2 T}{\partial x^2}+\dot{q},\quad 0<x<L.$$

- $k$ — conductivity, W/(m·K)
- $\rho$ — density, kg/m³
- $c_p$ — specific heat, J/(kg·K)
- $\dot{q}$ — internal volumetric heating, W/m³
- $T$ — temperature, °C (temperature *differences* also have units K)
- $\alpha=k/(\rho c_p)$ — thermal diffusivity, m²/s

All properties are fixed with respect to space and time; no phase changes are modeled.

## Cell-centered finite-volume discretization

Divide $L$ into $N$ equal cells with $\Delta x=L/N$. Integrate over cell $i$, per unit area normal to the slab:

$$\rho c_p\Delta x\,\frac{dT_i}{dt}=\frac{k}{\Delta x}(T_{i-1}-T_i)+\frac{k}{\Delta x}(T_{i+1}-T_i)+\dot{q}\Delta x.$$

At a boundary, substitute the corresponding outward heat flux, with *outward positive*. For example, at the **leftmost** cell:

$$\rho c_p\Delta x\,\frac{dT_0}{dt}=\frac{k}{\Delta x}(T_1-T_0)+\dot{q}\Delta x-q_{L,\mathrm{out}}.$$

Four types of boundaries are implemented:

| Type | Surface model | Outward heat flux |
|---|---|---|
| Insulated | $\partial T/\partial n=0$ | $q_{out}=0$ |
| Prescribed outward heat flux | $q_{out}=q_\mathrm{given}$ | $q_{out}=q_\mathrm{given}$ |
| Fixed surface temperature | Cell-to-surface half-cell resistance | $q_{out}=(2k/\Delta x)(T_{cell}-T_b)$ |
| Convection | Half-cell conduction + convection in series | $q_{out}=(T_{cell}-T_\infty)/[\Delta x/(2k)+1/h]$ |

This resistance treatment is important: simply applying $h(T_{cell}-T_\infty)$ does **not** account for the half-cell conduction between the cell center and the physical surface.

## Time integration

The spatial discretization yields a tridiagonal linear ODE system

$$\frac{d\mathbf{T}}{dt}=A\mathbf{T}+\mathbf{b}.$$

**Crank–Nicolson:**

$$(I-\tfrac{\Delta t}{2}A)\mathbf{T}^{n+1}=(I+\tfrac{\Delta t}{2}A)\mathbf{T}^{n}+\Delta t\,\mathbf{b}.$$

Implemented through SciPy's **banded tridiagonal solve** rather than a dense matrix inverse. For the linear problem it is A-stable, though overly large time steps can reduce accuracy and can produce nonphysical transient oscillations.

**Forward Euler:**

$$\mathbf{T}^{n+1}=\mathbf{T}^{n}+\Delta t(A\mathbf{T}^{n}+\mathbf{b}).$$

Explicit Euler is safeguarded with a conservative positivity condition: $\Delta t\leq 1/\max_i(-A_{ii})$. This is a sufficient stability/monotonicity limit for the present conduction discretization, not a claim about all PDE solvers.

## Discrete energy accounting

Stored energy relative to a fixed reference is evaluated per slab area:

$$E(t)=\rho c_p\Delta x\sum_iT_i.$$

Net incoming heat per unit area:

$$\dot Q_{net}=\dot q L-q_{L,out}-q_{R,out}.$$

For Crank–Nicolson, the heat input is accumulated using the trapezoidal rule to match the time integration. The reported residual is

$$R_E(t)=[E(t)-E(0)]-\int_0^t \dot Q_{net}(s)\,ds.$$

The residual should be near floating-point roundoff for correctly integrated simulations; this is a **discrete conservation check**, not empirical verification of material models.

## Analytical verification case

For an initially uniform slab with both surfaces held at the **same** temperature $T_b$, and no internal source, the Fourier series is

$$\frac{T(x,t)-T_b}{T_i-T_b}=\sum_{n=1,3,5,...}^{\infty}\frac{4}{n\pi}\sin\left(\frac{n\pi x}{L}\right)\exp\left(-\alpha\frac{n^2\pi^2}{L^2}t\right).$$

`thermalphysics verify` compares this analytical solution with the finite-volume solution on successively refined grids. The $t=0$ discontinuity at the physical boundary is treated separately.

## Explicit limitations

- One spatial dimension and uniform cell spacing
- One homogeneous layer, constant material properties
- No radiative boundary condition or nonlinear temperature-dependent properties
- No phase changes, ablation, pyrolysis, chemistry, fluid flow, or fluid–structure coupling
- Not a hypersonic environment model, heat-shield design tool, flight simulator, or certification software
- Does not model any real NASA/industry mission geometry, dataset or hardware

### Independent references

- Standard undergraduate/graduate derivation of the heat equation and Fourier-series solution; see a reputable heat-transfer textbook, e.g. *Fundamentals of Heat and Mass Transfer* (Incropera et al.).
- [SciPy documentation — `scipy.linalg.solve_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html).
