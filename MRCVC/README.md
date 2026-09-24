# MRCVC — Motor Rotativo de Combustión a Volumen Constante

This directory contains the computational fluid dynamics (CFD) configurations, mesh motion strategies, and numerical setups dedicated to simulating the **Motor Rotativo de Combustión a Volumen Constante (MRCVC)** [Constant-Volume Combustion Rotary Engine].

---

## Engine History, Patent & Academic Context

The **MRCVC** is a novel internal combustion engine design invented and patented in Argentina by researcher **Jorge A. Toth**. 

### Patent & Registration Details
* **Inventor:** Jorge A. Toth
* **Registry Body:** INPI (Instituto Nacional de la Propiedad Industrial), Argentina
* **Application / File No.:** `P 19960105411` (Original filing in 1996)
* **Patent Resolution:** `AR004806B1` (Definitively published in 2004)

### Research & Development
The modeling, design, and optimization of this rotary engine have been actively investigated within the **Department of Applied Mechanics** at the **Faculty of Engineering, Universidad Nacional del Comahue** (Neuquén, Argentina), supported by computational simulations and academic collaboration involving **CONICET** researchers. A group of engineering students and faculty advisors continues to build upon this work to evaluate its physical viability.

### Working Principle & Advantages
Unlike conventional reciprocating piston engines or the Wankel rotary engine, Toth's design achieves **constant-volume (isochoric) combustion**. Theoretically, this provides key thermodynamic and mechanical benefits:
* Higher thermal efficiency by reaching maximum thermodynamic availability.
* Reduced specific fuel consumption.
* Lower mechanical wear and reduced operational noise levels.

---

## Current Work & Methodology

Before coupling complex reacting flows or turbulence models, current numerical efforts focus on **verifying dynamic mesh stability** under large boundary displacements using OpenFOAM:

* **Parametric Mesh Generation:** The baseline 2D/3D geometry is constructed using `blockMeshDict`.
* **Motion Modeling:** Boundary node movements along the chamber walls (`inlet` and `outlet` geometric boundaries) are defined through custom C++ displacement logic in `0/pointDisplacement`.
* **Deformation Benchmarking:** Testing the operational limits of Laplacian motion solvers (`displacementLaplacian`) across varying rotor displacement angles to prevent cell skewness and mesh invalidation.

---

## 🐍 Python Prototyping (`MRCVCLaplaceMeshingDynamic.py`)

In addition to OpenFOAM cases, this directory includes `MRCVCLaplaceMeshingDynamic.py`, a standalone Python script used to prototype and visualize 2D structured mesh generation and deformation.

* **Methodology:** Discretizes Laplace's boundary value problem ($\nabla^2 x = 0, \nabla^2 y = 0$) using finite differences to smoothly distribute interior grid nodes based on moving rotor boundary profiles.
* **Purpose:** Serves as a fast numerical benchmark to verify boundary kinematics, node trajectories, and coordinate mapping before deploying full 3D OpenFOAM simulations.

## 📂 Minimal Case Structure (Mesh Motion Testing)

To keep the repository clean and lightweight during this testing phase, only the essentials required to execute `moveDynamicMesh` are included:

```text
MRCVC/
├── 0/
│   └── pointDisplacement    # Motion boundary condition definitions
├── constant/
│   └── dynamicMeshDict      # Motion solver specifications (Laplacian)
└── system/
    └── blockMeshDict        # Parametric mesh geometry definition
    ├── controlDict          # Time-step and execution parameters
    ├── fvSchemes            # Discretization schemes for point displacement
    └── fvSolution           # Linear matrix solvers for mesh movement
