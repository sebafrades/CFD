#include "Finite_Differences_Mesh.hpp"
#include "FV_Classes"

int main(){

    std::pair<Matrix,Matrix> mesh = generateMesh();
    Matrix x = mesh.first;
    Matrix y = mesh.second;

    // Create a vector of Points
    std::vector<Point> PointField;

    int rows = x.size();
    int cols = x[0].size();

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            Point p {x[i][j], y[i][j]}; // create object type P with (x,y)
            PointField.push_back(p);
        }
    }

    // Create a vector of faces
    std::vector<FV_face> Face_vector;

    for (int i = 0; i < rows ; ++i) {
        for(auto it = PointField.begin() + cols*i; it != PointField.begin() + cols*(i+1) - 1; ++it){
            //std::cout << "Point: (" << it->x << ", " << it->y << ")" << std::endl;
            FV_face face(*it, *(it + 1));
            Face_vector.push_back(face);
        }
    }

    for (int i = 0; i < rows - 1; ++i) {
        for(auto it = PointField.begin() + cols*i; it != PointField.begin() + cols*(i+1); ++it) {
            FV_face face(*it, *(it+cols));
            Face_vector.push_back(face);
        }
    }

    std::vector<FV_Cell> Cell_vector;

    for (auto it = Face_vector.begin(); it != Face_vector.end() - (cols-1)*rows - 1; ++it) {
        auto top = *it;
        auto right = *(it + (cols-1)*rows + 1);
        auto bottom = *(it + (cols-1));
        auto left = *(it + (cols-1)*rows);
        FV_Cell cell(top, right, bottom, left);
        Cell_vector.push_back(cell);
    }

}


