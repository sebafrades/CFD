# Elliptic Mesh Generation – Bump Case

This project implements a simple **elliptic grid generation method** in C++ and visualizes the resulting mesh using Python. The case consists of a 2D domain with a small geometric bump, serving as a basic exercise in structured mesh generation for CFD applications.

---

## Overview

* Solve elliptic mesh equations in **C++**
* Generate a smooth structured grid around a **bump geometry**
* Export mesh data to file
* Visualize the mesh using **Python**

The goal is to understand how elliptic PDE-based mesh generation improves grid smoothness compared to algebraic approaches.

---

## Method

The mesh is obtained by solving elliptic equations for grid coordinates, which iteratively smooth the mesh and distribute points more uniformly across the domain. This type of mesh generation is commonly used in CFD because mesh quality strongly affects solution accuracy and stability.

---

## Project Structure

```
Elliptic_meshes/
│
├── cpp/            # C++ solver (mesh generation)
├── python/         # Visualization scripts
├── data/           # Output mesh files
└── README.md
```

---

## How to Run

### 1. Compile the solver

```bash
g++ -O2 -o mesh main.cpp
```

### 2. Run

```bash
./mesh
```

This generates the mesh file in the output directory.

---

### 3. Plot the mesh

```bash
python plot.py
```

---

## Example Case

* Geometry: 2D domain with a **bump**
* Mesh type: structured grid
* Method: elliptic smoothing

The resulting mesh shows smoother spacing and better quality near the geometry compared to simple interpolation.

---

## Author

Sebastián Frades
Mechanical Engineering – CFD / Turbomachinery

