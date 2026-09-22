#ifndef FV_CELL_H
#define FV_CELL_H

#include<vector>

struct Point {
    double x;
    double y;
};

enum class FaceType { Interior, Inlet, Outlet, WallTop, WallBottom };

class FV_face {
public:
    FV_face(int p1_idx, int p2_idx)
        : p1_idx(p1_idx), p2_idx(p2_idx),
          type(FaceType::Interior),
          ownerCell(-1), neighborCell(-1),
          flux{0.0, 0.0, 0.0, 0.0} {}

    int getP1Index() const { return p1_idx; }
    int getP2Index() const { return p2_idx; }

    void setType(FaceType t) { type = t; }
    FaceType getType() const { return type; }

    void setOwner(int idx)    { ownerCell = idx; }
    void setNeighbor(int idx) { neighborCell = idx; }
    int getOwner() const    { return ownerCell; }
    int getNeighbor() const { return neighborCell; }

    void setFlux(const std::array<double,4>& f) { flux = f; }
    const std::array<double,4>& getFlux() const { return flux; }

private:
    int p1_idx;
    int p2_idx;
    FaceType type;
    int ownerCell;    // index into mesh.cells, -1 if unset
    int neighborCell; // index into mesh.cells, -1 if boundary face
    std::array<double,4> flux;
};

class FV_Cell {
public:
    FV_Cell(int ft, int fr, int fb, int fl)
        : faceTop(ft), faceRight(fr), faceBottom(fb), faceLeft(fl),
          rho(0.0), rhou(0.0), rhov(0.0), Et(0.0) {}

    int getFaceTop()    const { return faceTop; }
    int getFaceRight()  const { return faceRight; }
    int getFaceBottom() const { return faceBottom; }
    int getFaceLeft()   const { return faceLeft; }

    void setState(double rho_, double rhou_, double rhov_, double Et_) {
        rho = rho_; rhou = rhou_; rhov = rhov_; Et = Et_;
    }

    double getRho()  const { return rho; }
    double getRhoU() const { return rhou; }
    double getRhoV() const { return rhov; }
    double getEt()   const { return Et; }

    double getU() const { return rhou / rho; }
    double getV() const { return rhov / rho; }

    // --- geometry ---

    double area(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);

        std::vector<double> b = {p4.x - p1.x, p4.y - p1.y};
        std::vector<double> d = {p2.x - p3.x, p2.y - p3.y};
        return 0.5 * (b[0] * d[1] - b[1] * d[0]);
    }

    double delta_x_lef(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p1.x - p3.x;
    }

    double delta_x_rig(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p2.x - p4.x;
    }

    double delta_x_inf(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p4.x - p1.x;
    }

    double delta_x_sup(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p2.x - p3.x;
    }

    double delta_y_lef(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p1.y - p3.y;
    }

    double delta_y_rig(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p2.y - p4.y;
    }

    double delta_y_inf(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p4.y - p1.y;
    }

    double delta_y_sup(const std::vector<FV_face>& faces, const std::vector<Point>& points) const {
        Point p1, p2, p3, p4;
        getCorners(faces, points, p1, p2, p3, p4);
        return p2.y - p3.y;
    }

private:
    void getCorners(const std::vector<FV_face>& faces, const std::vector<Point>& points,
                     Point& p1, Point& p2, Point& p3, Point& p4) const {
        p1 = points[faces[faceTop].getP1Index()];
        p2 = points[faces[faceTop].getP2Index()];
        p3 = points[faces[faceBottom].getP1Index()];
        p4 = points[faces[faceBottom].getP2Index()];
    }

    int faceTop;
    int faceRight;
    int faceBottom;
    int faceLeft;

    double rho;
    double rhou;
    double rhov;
    double Et;
};

#endif // FV_CELL_H