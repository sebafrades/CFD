# MRCVC — Motor Rotativo de Combustión a Volumen Constante

This directory contains the computational fluid dynamics (CFD) configurations, mesh motion strategies, and numerical setups dedicated to simulating the **Motor Rotativo de Combustión a Volumen Constante (MRCVC)** [Constant-Volume Combustion Rotary Engine].

---

## 📌 Engine History & Context

The **MRCVC** is a non-conventional internal combustion engine conceived and patented by Argentine inventor **Eduardo Antonio Fernández**. The engine's core design aims to achieve constant-volume thermodynamic combustion within a rotary architecture, offering high thermal efficiency, reduced mechanical complexity, and continuous rotary motion.

Currently, a **collaborative team of engineering students and faculty advisors** (including work developed within the Department of Applied Mechanics) is working on the numerical modeling, dynamic mesh validation, and physical simulation of the MRCVC to evaluate its internal flow physics, rotor motion constraints, and thermodynamic behavior.

---

## 🔬 Current Work & Methodology

Before coupling complex reacting flows or turbulence models, current efforts focus on **verifying dynamic mesh stability** under large boundary displacements using OpenFOAM:

* **Parametric Mesh Generation:** The baseline 2D/3D geometry is constructed using `blockMeshDict`.
* **Motion Modeling:** Boundary node movements along the chamber walls (`inlet` and `outlet` geometric boundaries) are defined through custom C++ displacement logic in `0/pointDisplacement`.
* **Deformation Benchmarking:** Testing the operational limits of Laplacian motion solvers (`displacementLaplacian`) across varying rotor displacement angles to prevent cell skewness and mesh invalidation.

---

## 📂 Minimal Case Structure (Mesh Motion Testing)

To keep the repository clean and efficient during this testing phase, only the essentials for running `moveDynamicMesh` are included:

```text
MRCVC/
├── 0/
│   └── pointDisplacement    # Motion boundary condition definitions
├── constant/
│   ├── dynamicMeshDict      # Motion solver specifications (Laplacian)
│   └── polyMesh/
│       └── blockMeshDict    # Parametric mesh geometry definition
└── system/
    ├── controlDict          # Time-step and execution parameters
    ├── fvSchemes            # Discretization schemes for point displacement
    └── fvSolution           # Linear matrix solvers for mesh movement
