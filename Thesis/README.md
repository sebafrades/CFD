# Numerical Resolution of Flow Through a Sudden Symmetric Expansion

This repository contains the code, simulated case configurations, and the complete thesis document for the Professional Integrative Project submitted for the degree of **Mechanical Engineer** at the National University of Comahue (Neuquén, Argentina, 2024).

* **Author:** Sebastián Alexis Frades
* **Advisor:** Dr. Ricardo Adolfo Prado
* **Collaborator:** Dr. Ezequiel López

---

## 📌 Project Overview

The primary objective of this work is the numerical resolution of the governing equations of fluid mechanics in conduits featuring abrupt geometric changes (sudden symmetric expansion). The complexity of the problem is approached incrementally, starting from two-dimensional potential and laminar flows to three-dimensional turbulent cases and sediment transport (dispersed/particulate phase).

The project emphasizes **verification** and **validation** of computational modeling by comparing custom-built finite-difference codes with **OpenFOAM** (finite volume method) simulations against analytical, numerical, and experimental data from literature.

---

## 🛠️ Methodology & Case Studies

1. **Two-Dimensional Cases (Finite Differences & OpenFOAM):**
   * 2D Potential flow (Uniform and parabolic inlet profiles).
   * 2D Viscous laminar flow (Stream Function - Vorticity $\psi-\zeta$ formulation, using Centered and 3rd-order UPWIND schemes)[cite: 1].
   * Multigrid acceleration for higher Reynolds number computations[cite: 1].
   * Axisymmetric 2D viscous laminar flow[cite: 1].
   * Validation of the laminar solver in OpenFOAM (`icoFoam`)[cite: 1].
2. **Three-Dimensional Turbulent Cases (OpenFOAM):**
   * RANS modeling of axisymmetric turbulent flows[cite: 1].
   * Qualitative and quantitative comparison of two-equation turbulence models ($k-\epsilon$ and $k-\omega$) against experimental data (Khezzar et al., 1986)[cite: 1].
3. **Sediment Transport (OpenFOAM):**
   * Implementation of an Eulerian-Lagrangian approach to evaluate particle dynamics and critical erosion zones in sudden expansions[cite: 1].

---

## 📂 Repository Structure

```text
.
├── 2D_FiniteDifferences/   # Custom C++/Python/MATLAB codes for finite differences (potential, psi-zeta, multigrid)
├── OpenFOAM/
│   ├── Laminar_2D/        # OpenFOAM cases using icoFoam (benchmarked against finite differences)
│   ├── Turbulent_3D/      # RANS cases using k-epsilon and k-omega models
│   └── Sediment_Transport/# Multiphase flow / Lagrangian particle tracking configurations
├── Thesis_Document.pdf    # Full thesis PDF document
└── README.md              # General repository documentation
