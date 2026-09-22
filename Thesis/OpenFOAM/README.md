# OpenFOAM Simulation Cases

This directory contains cases discretized via the **finite volume method** using the open-source platform **OpenFOAM**[cite: 1].

## 📂 Directory Layout

```text
OpenFOAM/
├── Laminar_2D/          # 2D laminar viscous flow computed with icoFoam for direct comparison with FD
├── Turbulent_3D/        # 3D axisymmetric mesh evaluated with simpleFoam
│   ├── kEpsilon/        # Configuration and results for the k-epsilon model
│   └── kOmega/          # Configuration and results for the k-omega model
└── Sediment_Transport/  # Lagrangian sediment transport over background turbulent flow
