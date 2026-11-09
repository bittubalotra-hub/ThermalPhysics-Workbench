"""Optional interactive demonstrator: run `streamlit run app.py`.

All inputs are illustrative and do not describe real flight hardware.
"""
import matplotlib.pyplot as plt
import streamlit as st
from thermalphysics import Boundary, Material, Slab, SolverOptions, simulate

st.set_page_config(page_title="ThermalPhysics Workbench", page_icon="🌡️", layout="wide")
st.title("ThermalPhysics Workbench")
st.caption("Independent research demo · 1D transient conduction · generic materials · SI units")

with st.sidebar:
    st.header("Slab & material")
    thickness_mm = st.slider("Thickness (mm)", 5, 60, 20, step=1)
    cells = st.slider("Finite-volume cells", 20, 200, 80, step=10)
    k = st.slider("Conductivity (W/m/K)", 1.0, 150.0, 45.0)
    density = st.number_input("Density (kg/m³)", 500.0, 20000.0, 7800.0, 100.0)
    heat_capacity = st.number_input("Heat capacity (J/kg/K)", 100.0, 3000.0, 500.0, 50.0)
    st.header("Boundary conditions")
    initial = st.slider("Initial temperature (°C)", 0.0, 300.0, 180.0)
    ambient = st.slider("Ambient temperature (°C)", -20.0, 100.0, 22.0)
    left_h = st.slider("Left convection h (W/m²/K)", 5.0, 500.0, 55.0)
    right_h = st.slider("Right convection h (W/m²/K)", 5.0, 500.0, 120.0)
    st.header("Integration")
    duration = st.slider("Total duration (s)", 20, 1200, 450, step=10)
    timestep = st.selectbox("Time step (s)", [0.1, 0.25, 0.5, 1.0], index=2)

material = Material(k, density, heat_capacity, "User-defined illustrative solid")
model = Slab(thickness_mm / 1000, cells, material,
             Boundary.convection(left_h, ambient), Boundary.convection(right_h, ambient), initial)
options = SolverOptions(timestep, int(duration/timestep), max(1, int(duration/timestep)//35))
result = simulate(model, options)

c1, c2, c3 = st.columns(3)
c1.metric("Thermal diffusivity", f"{material.diffusivity_m2_s:.3e} m²/s")
c2.metric("Final mean temperature", f"{result.final_temperature_c.mean():.2f} °C")
c3.metric("Energy balance residual", f"{result.max_abs_balance_error_j_m2:.2e} J/m²")

fig, ax = plt.subplots(figsize=(10, 5))
for idx in [0, len(result.times_s)//4, len(result.times_s)//2, -1]:
    ax.plot(result.positions_m*1e3, result.temperatures_c[idx], lw=2,
            label=f"t={result.times_s[idx]:g} s")
ax.set(xlabel="Depth (mm)", ylabel="Temperature (°C)", title="Temperature through the slab")
ax.grid(alpha=0.25)
ax.legend()
st.pyplot(fig, clear_figure=True)

st.info("Educational model only: 1D, homogeneous material, constant thermal properties, no radiation or phase change. Not suitable for flight-hardware qualification.")
