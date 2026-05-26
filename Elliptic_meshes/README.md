Structured mesh generation for 2D domains, covering algebraic and elliptic methods, iterative solvers, and convergence analysis.

TP6 — Algebraic and Elliptic Meshes in Channel Domains
Applies algebraic mesh generation (trivial transformation ξ=x, η=y\xi = x,\ \eta = y
ξ=x, η=y) to two channel geometries: one with a rectangular step and one with a parabolic upper boundary. The effect of node spacing (Δx=Δy=a\Delta x = \Delta y = a
Δx=Δy=a and a/2a/2
a/2) on the resulting mesh is analyzed.
A second problem introduces elliptic mesh generation via the Poisson equation on a trapezoidal domain, where source terms control mesh clustering. Two different initial meshes are used to study convergence, showing mesh redistribution and log⁡∣E∣\log|E|
log∣E∣ vs. iteration plots.
TP7 — Elliptic Mesh Generation for a Duct with a Circular Arc Bump
Generates meshes for a 2D duct with a circular arc protrusion using both Laplace (P=Q=0P=Q=0
P=Q=0) and Poisson (P,Q≠0P,Q \neq 0
P,Q=0) elliptic equations. Two boundary node distributions are considered: uniform spacing (49×17 grid) and algebraic wall-normal clustering (65×17 grid).
The algebraic systems are solved with Gauss-Seidel, successive line relaxation by rows, and successive line relaxation by columns. The optimal over-relaxation factor ω\omega
ω is determined for each case along with the iteration count to convergence.
