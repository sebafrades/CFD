#include "mesh.hpp"

int main(){

    std::pair<Matrix,Matrix> mesh = generateMesh();
    Matrix x = mesh.first;
    Matrix y = mesh.second;

    saveMatrix(x, y, "grid.csv");
}