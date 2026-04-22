#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

using Matrix = vector<vector<double>>;

Matrix createMatrix(int rows, int cols, char type);
void iterateMatrix(Matrix& x, Matrix& y, double error);
double algebraicFunction(double i);
void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename);

constexpr int N = 65; // xi
constexpr int M = 17; // eta

int main(){

    Matrix x = createMatrix(M, N, 'x');
    Matrix y = createMatrix(M, N, 'y');

    iterateMatrix(x, y, 1e-6);

    saveMatrix(x, y, "grid.csv");
}

// -----------------------------
// Algebraic stretching function
// -----------------------------
double algebraicFunction(double i) {
    double K = 1.0864;
    int L = 17;
    return (pow(K,i-1) - 1.0)/(pow(K,L-1) - 1.0);
}

// -----------------------------
// Boundary generation
// -----------------------------
Matrix createMatrix(int rows, int cols, char type) {

    Matrix matrix(rows, vector<double>(cols, 0.0));

    if(type == 'x'){

        for(int j = 0; j < cols; j++){

            double x_val;

            if (j < 17){
                int idx = j + 1;
                double s = algebraicFunction(17 - idx + 1);
                x_val = 1.0 * (1.0 - s);
            }
            else if (j < 49){
                x_val = 1.0 + (double(j - 17)/32.0);
            }
            else{
                int idx = j + 1;
                double s = algebraicFunction(idx - 48);
                x_val = 2.0 + s;
            }

            matrix[0][j] = x_val;
            matrix[rows-1][j] = x_val;
        }

        // vertical boundaries (straight)
        for(int i = 0; i < rows; i++){
            matrix[i][0] = 0.0;
            matrix[i][cols-1] = 3.0;
        }
    }

    // -----------------------------
    if(type == 'y'){

        // TOP boundary
        for(int j = 0; j < cols; j++){
            matrix[0][j] = 1.0;
        }

        // BOTTOM boundary (parabolic bump)
        for(int j = 0; j < cols; j++){

            if (j < 17){
                matrix[rows-1][j] = 0.0;
            }
            else if (j < 49){

                // reconstruct x consistently
                double x_val = 1.0 + (double(j - 17)/32.0);

                double xc = 1.5;
                double half_width = 0.5;
                double H = 0.2;

                double xi = (x_val - xc)/half_width;

                matrix[rows-1][j] = H * (1.0 - xi*xi);
            }
            else{
                matrix[rows-1][j] = 0.0;
            }
        }

        // SIDE boundaries (stretched)
        for(int i = 0; i < rows; i++){
            double val = algebraicFunction(i+1);

            matrix[rows-1-i][0] = val;
            matrix[rows-1-i][cols-1] = val;
        }
    }

    return matrix;
}

// -----------------------------
// Elliptic solver with P,Q
// -----------------------------
void iterateMatrix(Matrix& x, Matrix& y, double error) {

    int rows = x.size();
    int cols = x[0].size();

    Matrix Phi(rows, vector<double>(cols, 0.0));
    Matrix Psi(rows, vector<double>(cols, 0.0));

    double residual = 1.0;

    const double w = 1;     // relaxation
    const double damping = 1; // P,Q damping

    while (residual > error){

        residual = 0.0;

        // ---- Phi (top & bottom)
        for (int j = 1; j < cols-1; j++){
            for (int i : {0, rows-1}){

                double x_xi = (x[i][j+1] - x[i][j-1]) / 2.0;
                double y_xi = (y[i][j+1] - y[i][j-1]) / 2.0;

                double x_xixi = x[i][j+1] - 2*x[i][j] + x[i][j-1];
                double y_xixi = y[i][j+1] - 2*y[i][j] + y[i][j-1];

                double denom = x_xi*x_xi + y_xi*y_xi + 1e-12;

                Phi[i][j] = -(x_xixi*x_xi + y_xixi*y_xi) / denom;
            }
        }

        // ---- Psi (left & right)
        for (int i = 1; i < rows-1; i++){
            for (int j : {0, cols-1}){

                double x_eta = (x[i+1][j] - x[i-1][j]) / 2.0;
                double y_eta = (y[i+1][j] - y[i-1][j]) / 2.0;

                double x_etaeta = x[i+1][j] - 2*x[i][j] + x[i-1][j];
                double y_etaeta = y[i+1][j] - 2*y[i][j] + y[i-1][j];

                double denom = x_eta*x_eta + y_eta*y_eta + 1e-12;

                Psi[i][j] = -(x_etaeta*x_eta + y_etaeta*y_eta) / denom;
            }
        }

        // ---- Interpolation
        for (int i = 1; i < rows-1; i++){
            for (int j = 1; j < cols-1; j++){

                Phi[i][j] = (1.0 - double(i)/(rows-1)) * Phi[0][j]
                          + (double(i)/(rows-1)) * Phi[rows-1][j];

                Psi[i][j] = (1.0 - double(j)/(cols-1)) * Psi[i][0]
                          + (double(j)/(cols-1)) * Psi[i][cols-1];
            }
        }

        // ---- Main iteration
        for (int i = 1; i < rows-1; i++){
            for (int j = 1; j < cols-1; j++){

                double x_eta = (x[i+1][j] - x[i-1][j]) / 2.0;
                double x_xi  = (x[i][j+1] - x[i][j-1]) / 2.0;
                double y_eta = (y[i+1][j] - y[i-1][j]) / 2.0;
                double y_xi  = (y[i][j+1] - y[i][j-1]) / 2.0;

                double x_xi_eta = 0.25*(x[i+1][j+1] - x[i+1][j-1]
                                     - x[i-1][j+1] + x[i-1][j-1]);

                double y_xi_eta = 0.25*(y[i+1][j+1] - y[i+1][j-1]
                                     - y[i-1][j+1] + y[i-1][j-1]);

                double alpha = x_eta*x_eta + y_eta*y_eta;
                double beta  = x_xi*x_eta + y_xi*y_eta;
                double gamma = x_xi*x_xi + y_xi*y_xi;

                double J = x_xi*y_eta - x_eta*y_xi;

                double P = damping * Phi[i][j] * (x_xi*x_xi + y_xi*y_xi);
                double Q = damping * Psi[i][j] * (x_eta*x_eta + y_eta*y_eta);

                double x_new = (
                    2*beta*x_xi_eta
                    - alpha*(x[i][j+1] + x[i][j-1])
                    - gamma*(x[i+1][j] + x[i-1][j])
                    - J*J*(P*x_xi + Q*x_eta)
                ) / (-2*(alpha + gamma));

                double y_new = (
                    2*beta*y_xi_eta
                    - alpha*(y[i][j+1] + y[i][j-1])
                    - gamma*(y[i+1][j] + y[i-1][j])
                    - J*J*(P*y_xi + Q*y_eta)
                ) / (-2*(alpha + gamma));

                residual = max(residual, fabs(x_new - x[i][j]));

                // relaxation
                x[i][j] = (1-w)*x[i][j] + w*x_new;
                y[i][j] = (1-w)*y[i][j] + w*y_new;
            }
        }

        cout << "Residual: " << residual << endl;
    }
}

// -----------------------------
void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename) {
    ofstream file(filename);

    for (int i = 0; i < x.size(); i++){
        for (int j = 0; j < x[0].size(); j++){
            file << x[i][j] << "," << y[i][j];
            if (j < x[0].size()-1) file << ",";
        }
        file << "\n";
    }
}