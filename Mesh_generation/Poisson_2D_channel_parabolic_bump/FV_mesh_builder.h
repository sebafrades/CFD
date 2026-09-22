#ifndef FV_MESH_BUILDER_H
#define FV_MESH_BUILDER_H

#include <vector>
#include "Finite_Differences_Mesh.hpp"
#include "FV_Classes.h"

struct FVMesh {
    std::vector<Point> points;
    std::vector<FV_face> faces;
    std::vector<FV_Cell> cells;
};

FVMesh buildFVMesh(const Matrix& x, const Matrix& y);

#endif // FV_MESH_BUILDER_H