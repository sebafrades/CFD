#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

using Matrix = vector<vector<double>>;

// matrix creation function declaration
Matrix createMatrix(int rows, int cols, char type);

// matrix iteration process function declaration
void iterateMatrix(Matrix& x, Matrix& y, double error);

// algebraic function
double algebraicFunction(double i);

constexpr int N = 65;
constexpr int M = 17;

void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename);

int main(){

    Matrix x = createMatrix(M, N, 'x');
    Matrix y = createMatrix(M, N, 'y');

    //for (int j = 0; j < N; j++) {
    //    cout << x[M-1][j] << " ";
    //}
    //cout << endl;

    iterateMatrix(x, y, 1e-7);

    saveMatrix(x, y, "grid.csv");
}

// function definition
Matrix createMatrix(int rows, int cols, char type) {
    
    Matrix matrix(rows, vector<double>(cols, 0.0));

    switch(type){
        case 'x': 

            // superior and inferior row boundary
            for(int i = 0; i < cols; i++){

                if (i <= 16){  // first 17 points (0–16 or 1-17 in 1-based indexing)

                    int idx = i + 1;  // maps: i=0→1, i=16→17
                    double s = algebraicFunction(17 - idx + 1);

                    matrix[0][i] = 1.0 * (1.0 - s);
                    matrix[rows-1][i] = matrix[0][i];
                }

                else if (i >= 16 && i <= 48){  // middle region (17–48)

                    matrix[0][i] = 1.0 + (1.0)*(double(i - 17)/32.0);
                    matrix[rows-1][i] = matrix[0][i];
                }

                else if (i >= 48 && i <= N-1){  // last 17 points (49–64)

                    int idx = i + 1;  // maps: i=48→49, i=64→65
                    double s = algebraicFunction(idx - 48);

                    matrix[0][i] = (1.0 + 0.5 + 0.5) + 1.0 * s;
                    matrix[rows-1][i] = matrix[0][i];
                }
            }

            // left and right column boundary
            for(int i=0; i<rows; i++){
                matrix[i][0] = 0.0;
                matrix[i][cols-1] = matrix[rows-1][N-1];
            }

            break;

        case 'y':
            // y boundary
            for(int i = 0; i < cols; i++){
                matrix[0][i] = 1.0;

                if (i <= 16) {
                    matrix[rows-1][i] = 0.0;
                }
                else if (i >= 16 && i <= 48){

                    double x_val = 1.0 + (1.0)*(double(i - 17)/32.0);

                    double xc = 1.5;
                    double half_width = 0.5;
                    double H = 0.2;

                    double xi = (x_val - xc) / half_width;

                    matrix[rows-1][i] = H * (1.0 - xi*xi);
                }
                else if (i >= 48 && i <= N-1){
                    matrix[rows-1][i] = 0.0;
                }
            }

            // columna izquierda y derecha
            for(int i=0; i<rows; i++){
                matrix[rows-1-i][0] = 1.0*algebraicFunction(i+1);
                matrix[rows-1-i][cols-1] = 1.0*algebraicFunction(i+1);
            } 
            break;
    }

    return matrix;
}

void iterateMatrix(Matrix& x, Matrix& y, double error) {

    double residual = 1.0;

    int rows = x.size();   // j (eta)
    int cols = x[0].size(); // i (xi)

    // Phi and Psi storage
    Matrix Phi(rows, vector<double>(cols, 0.0));
    Matrix Psi(rows, vector<double>(cols, 0.0));

    while (residual > error){

        residual = 0.0;

        // --- STEP 1: Compute Phi on top/bottom boundaries ---
        for (int j = 1; j < cols-1; j++){
            for (int i : {0, rows-1}){

                double x_xi = (x[i][j+1] - x[i][j-1]) / 2.0;
                double y_xi = (y[i][j+1] - y[i][j-1]) / 2.0;

                double x_xixi = x[i][j+1] - 2*x[i][j] + x[i][j-1];
                double y_xixi = y[i][j+1] - 2*y[i][j] + y[i][j-1];

                double denom = pow(x_xi, 2) + pow(y_xi, 2) + 1e-12;

                Phi[i][j] = -(x_xixi*x_xi + y_xixi*y_xi) / denom;
            }
        }

        // --- STEP 2: Compute Psi on left/right boundaries ---
        for (int i = 1; i < rows-1; i++){
            for (int j : {0, cols-1}){

                double x_eta = (x[i+1][j] - x[i-1][j]) / 2.0;
                double y_eta = (y[i+1][j] - y[i-1][j]) / 2.0;

                double x_etaeta = x[i+1][j] - 2*x[i][j] + x[i-1][j];
                double y_etaeta = y[i+1][j] - 2*y[i][j] + y[i-1][j];

                double denom = pow(x_eta, 2) + pow(y_eta, 2) + 1e-12;

                Psi[i][j] = -(x_etaeta*x_eta + y_etaeta*y_eta) / denom;
            }
        }

        // --- STEP 3: Interpolate Phi and Psi into interior ---
        for (int i = 1; i < rows-1; i++){
            for (int j = 1; j < cols-1; j++){

                // linear interpolation
                Phi[i][j] = (1.0 - double(i)/(rows-1)) * Phi[0][j]
                          + (double(i)/(rows-1)) * Phi[rows-1][j];

                Psi[i][j] = (1.0 - double(j)/(cols-1)) * Psi[i][0]
                          + (double(j)/(cols-1)) * Psi[i][cols-1];
            }
        }

        // --- STEP 4: Main elliptic iteration ---
        for (int i = 1; i < rows-1; i++){
            for (int j = 1; j < cols-1; j++){
                
                double x_eta = (x[i+1][j] - x[i-1][j]) / 2.0;
                double x_xi  = (x[i][j+1] - x[i][j-1]) / 2.0;
                double y_eta = (y[i+1][j] - y[i-1][j]) / 2.0;
                double y_xi  = (y[i][j+1] - y[i][j-1]) / 2.0;

                double x_xi_eta = 0.25 * (x[i+1][j+1] - x[i+1][j-1]
                                       - x[i-1][j+1] + x[i-1][j-1]);

                double y_xi_eta = 0.25 * (y[i+1][j+1] - y[i+1][j-1]
                                       - y[i-1][j+1] + y[i-1][j-1]);

                double alpha = pow(x_eta, 2) + pow(y_eta, 2);
                double beta  = x_xi*x_eta + y_xi*y_eta;
                double gamma = pow(x_xi, 2) + pow(y_xi, 2);

                double J = x_xi*y_eta - x_eta*y_xi;

                double P = Phi[i][j] * (x_xi*x_xi + y_xi*y_xi);
                double Q = Psi[i][j] * (x_eta*x_eta + y_eta*y_eta);

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

                double diff = fabs(x_new - x[i][j]);
                if (diff > residual) residual = diff;

                x[i][j] = x_new;
                y[i][j] = y_new;
            }
        }
    }
}

double algebraicFunction(double i) {
    double K = 1.0864;
    int L = 17;

    double f = (pow(K,i-1) - 1)/(pow(K,L-1) - 1);

    return f;
}

void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename) {
    std::ofstream file(filename);

    int rows = x.size();
    int cols = x[0].size();

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            file << x[i][j] << "," << y[i][j];
            if (j < cols-1) file << ",";
        }
        file << "\n";
    }
}
