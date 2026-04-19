# CFD

This repository contains a set of CFD exercises and implementations written mainly in C++ and Python, focused on structured grids and finite difference methods. The idea is to build everything from scratch to better understand the numerical and physical aspects behind the methods.

## What’s in here

* Elliptic mesh generation (Laplace / Poisson)
* Classical unsteady problems:

  * Second Stokes problem
  * Transient Poiseuille flow
  * Transient Burgers equation


All solvers are written in C++, and results are post-processed with Python.

## Approach

The focus is on:

* writing the numerical methods directly
* keeping the code simple and readable
* understanding convergence, stability, and boundary conditions

This is not meant to be production code — it’s closer to a numerical lab notebook.

## Author

Sebastián Frades
Mechanical Engineering
