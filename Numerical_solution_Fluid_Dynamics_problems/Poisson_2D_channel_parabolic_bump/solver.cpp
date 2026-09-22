#include "Finite_Differences_Mesh.hpp"
#include "FV_mesh_builder.h"
#include <iostream>

// ---------------------------------------------------------
// Freestream / reference conditions
// ---------------------------------------------------------
// Freestream conditions (nondimensional)
constexpr double gamma = 1.4;
constexpr double R     = 1.0;
constexpr double Minf  = 0.5;
constexpr double rhoInf = 1.0;
constexpr double pInf   = 1.0;
constexpr double Tinf   = pInf / (rhoInf * R);       // = 1.0
constexpr double aInf   = std::sqrt(gamma * R * Tinf);
constexpr double uInf   = Minf * aInf;
constexpr double vInf   = 0.0;

constexpr double Cv     = R / (gamma - 1.0);
constexpr double eInf   = Cv * Tinf;                                  // internal energy per unit mass
constexpr double EtInf  = rhoInf * (eInf + 0.5*(uInf*uInf + vInf*vInf)); // total energy per unit volume

// ---------------------------------------------------------
// Initialization
// ---------------------------------------------------------
void initializeField(FVMesh& mesh) {
    for (auto it = mesh.cells.begin(); it != mesh.cells.end(); ++it) {
        it->setState(rhoInf, rhoInf*uInf, rhoInf*vInf, EtInf);
    }
}
// Sets every cell's conservative state to freestream values

// ---------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------
void applyBoundaryConditions(FVMesh& mesh) {
    for (auto& face : mesh.faces) {
        switch (face.getType()) {
            case FaceType::WallTop:
            case FaceType::WallBottom:
                applyWallBC(face, /* owning cell */);
                break;
            case FaceType::Inlet:
                applyInletBC(face, /* owning cell */);
                break;
            case FaceType::Outlet:
                applyOutletBC(face, /* owning cell */);
                break;
            default:
                break; // interior face, nothing to do here
        }
    }
}

// ---------------------------------------------------------
// Flux computation
// ---------------------------------------------------------
// Given two neighboring cells' states (or a cell + boundary/ghost state),
// compute the averaged interface state, then F/G, then the combined
// normal flux for a face.
void computeFaceFluxes(FVMesh& mesh);

// ---------------------------------------------------------
// Residual assembly
// ---------------------------------------------------------
// For each cell, sum contributions from its four faces
// (per the F*dy - G*dx style equation)
void computeResidual(FVMesh& mesh /*, whatever residual storage you choose */);

// ---------------------------------------------------------
// Time marching
// ---------------------------------------------------------
void timeStep(FVMesh& mesh, double dt);
// Advances each cell's conservative state using its residual

// ---------------------------------------------------------
// Output
// ---------------------------------------------------------
void saveResults(const FVMesh& mesh, const std::string& filename);
// Write rho, u, v, p, M (or whatever you want) per cell to CSV

// ---------------------------------------------------------
int main(){

    // 1. Get the Poisson grid
    std::pair<Matrix,Matrix> gridResult = generateMesh();

    // 2. Build the FV mesh (points, faces, cells)
    FVMesh mesh = buildFVMesh(gridResult.first, gridResult.second);

    // 3. Initialize the field
    initializeField(mesh, fs);

    // 4. Main iteration loop
    double residualNorm = 1.0;
    double error = 1e-6;
    int iteration = 0;

    while (residualNorm > error){

        // apply boundary conditions
        // compute face fluxes
        // compute residual per cell
        // advance in time
        // update residualNorm
        // (print progress)

        iteration++;
    }

    // 6. Save results
    saveResults(mesh, "solution.csv");

    std::cout << "Done. Iterations: " << iteration << std::endl;
}