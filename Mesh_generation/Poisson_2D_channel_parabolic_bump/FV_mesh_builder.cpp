#include "FV_mesh_builder.h"

FVMesh buildFVMesh(const Matrix& x, const Matrix& y) {

    FVMesh mesh;

    int rows = x.size();
    int cols = x[0].size();

    // Points
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            Point p {x[i][j], y[i][j]};
            mesh.points.push_back(p);
        }
    }

    // Horizontal faces (within each row): connects (i,j) to (i,j+1)
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols - 1; ++j) {
            int p1_idx = i*cols + j;
            int p2_idx = i*cols + j + 1;

            FV_face face(p1_idx, p2_idx);

            if (i == 0) {
                face.setType(FaceType::WallTop);
            } else if (i == rows - 1) {
                face.setType(FaceType::WallBottom);
            }
            // else: leave as default Interior

            mesh.faces.push_back(face);
        }
    }

    // Vertical faces (between rows): connects (i,j) to (i+1,j)
    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 0; j < cols; ++j) {
            int p1_idx = i*cols + j;
            int p2_idx = (i+1)*cols + j;

            FV_face face(p1_idx, p2_idx);

            if (j == 0) {
                face.setType(FaceType::Inlet);
            } else if (j == cols - 1) {
                face.setType(FaceType::Outlet);
            }
            // else: leave as default Interior

            mesh.faces.push_back(face);
        }
    }

    // Cells
    int horizontalFaceCount = rows * (cols - 1);

    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 0; j < cols - 1; ++j) {

            int topIdx    = i*(cols-1) + j;
            int bottomIdx = (i+1)*(cols-1) + j;
            int leftIdx   = horizontalFaceCount + i*cols + j;
            int rightIdx  = horizontalFaceCount + i*cols + (j+1);

            FV_Cell cell(topIdx, rightIdx, bottomIdx, leftIdx);
            int cellIdx = static_cast<int>(mesh.cells.size());
            mesh.cells.push_back(cell);

            // assign owner/neighbor on each of this cell's faces
            for (int faceIdx : {topIdx, rightIdx, bottomIdx, leftIdx}) {
                if (mesh.faces[faceIdx].getOwner() == -1) {
                    mesh.faces[faceIdx].setOwner(cellIdx);
                } else {
                    mesh.faces[faceIdx].setNeighbor(cellIdx);
                }
            }
        }
    }

    return mesh;
}
