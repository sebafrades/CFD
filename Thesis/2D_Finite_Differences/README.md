# 2D Finite Differences Implementation (`2D_Finite_Differences`)

This directory contains custom numerical implementations using the **Finite Difference Method (FDM)** to analyze fluid flow through a sudden symmetric expansion[cite: 1].

## 📂 Directory Layout & Modules

* **`2D_Inviscid/`**: Solves two-dimensional irrotational/potential flow ($\psi$) using centered differences with uniform boundary profiles[cite: 1].
* **`2D_Inviscid_I_O_viscous/`**: Potential flow formulation incorporating viscous inlet and outlet profiles[cite: 1].
* **`2D_Viscous/`**: Solves the full 2D incompressible Navier-Stokes equations using the Stream Function - Vorticity ($\psi-\zeta$) formulation (including centered and UPWIND 3rd-order discretizations, as well as Multigrid acceleration)[cite: 1].
* **`Axisymmetric/`**: Implementation of the Stokes stream function equations for 2D axisymmetric expansion geometries in cylindrical coordinates[cite: 1].

## Usage & Execution

1. Navigate to the desired folder (e.g., `cd 2D_Viscous`).
2. Set domain dimensions ($H/h$, $L_1/h$, $L_2/h$) and operating parameters such as Reynolds number ($Re$)[cite: 1].
3. Run the primary solver scripts to compute stream functions, vorticity fields, pressure distribution ($C_p$), and reattachment length ($x_r$)[cite: 1].
