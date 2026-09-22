#ifndef MESH_HPP
#define MESH_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <utility>

using Matrix = std::vector<std::vector<double>>;

Matrix createMatrix(int rows, int cols, char type);
void iterateMatrix(Matrix& x, Matrix& y, double error);
double algebraicFunction(double i);
void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename);
std::pair<Matrix,Matrix> generateMesh();

#endif // MESH_HPP