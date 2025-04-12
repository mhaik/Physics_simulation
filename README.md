# Physics_simulation

## 🧊 Task 1 – Heat transfer in a hollow pipe

This project simulates the **steady-state temperature distribution** and **heat fluxes** in an infinitely long square pipe with an off-center rectangular hollow section. The pipe is surrounded by an external medium with different constant temperatures on each side (`T1`, `T2`, `T3`, `T4`) and exchanges heat with a hot internal fluid at temperature `T₀`. 

The simulation solves the **2D heat conduction equation** using **second-kind (Robin) boundary conditions** with known heat transfer coefficients (`λi` for the inside, `λ₀` for the outside). The results include:
- The **temperature field** `T(x, y)`
- The **heat flux vectors**
- **Verification of solution accuracy**, e.g. via grid refinement

---

## 🪐 Task 2 – Orbital simulation of planetoids

This simulation models the **orbital motion of planetoids** between the orbits of Mars and Jupiter. Assuming all orbits are **circular** and **coplanar**, and that the planetoids do not influence the motion of the planets, the system simplifies to a two-body (Sun-planetoid) problem.

The orbits are simulated for the following radii (in AU):
- `r₁ = 3.000 AU`
- `r₂ = 3.276 AU`
- `r₃ = 3.700 AU`

Using **Kepler’s laws**, the corresponding orbital velocities are computed, and the trajectories are visualized in a **heliocentric coordinate system**, giving a simplified model of motion in the asteroid belt region.
