---

### 2. `README.md` for `2D_FiniteDifferences/` (Custom Finite Difference Codes)

```markdown
# 2D Finite Difference Modules (Custom Codes)

This directory contains custom implementations based on the **finite difference method** to solve two-dimensional flows through sudden symmetric expansions[cite: 1].

## 📑 Directory Contents

* `potential_flow/`: Irrotational/inviscid flow calculations (Laplace equation for $\psi$) using $O(\Delta x^2)$ centered differences with uniform and parabolic inlet velocity profiles[cite: 1].
* `stream_vorticity_centered/`: 2D incompressible Navier-Stokes solver using the Stream Function - Vorticity ($\psi-\zeta$) formulation with centered schemes[cite: 1].
* `stream_vorticity_upwind/`: Implementation of 3rd-order Upwind scheme for the discretization of convective terms in the vorticity equation[cite: 1].
* `multigrid_acceleration/`: Multigrid acceleration algorithm (V-Cycle) designed to smooth high-frequency errors and accelerate convergence at higher Reynolds numbers (up to $Re = 360$)[cite: 1].
* `axisymmetric_laminar/`: Stokes stream function formulation for laminar flows in cylindrical/axisymmetric coordinates[cite: 1].

## 🚀 Getting Started

1. Select the subdirectory for the target case.
2. Define geometric parameters ($H/h$, $L_1/h$, $L_2/h$) and operating conditions (Reynolds number $Re$)[cite: 1].
3. Run the simulation scripts and use post-processing scripts to output streamlines, pressure coefficients $C_p$, and reattachment lengths $x_r$[cite: 1].
