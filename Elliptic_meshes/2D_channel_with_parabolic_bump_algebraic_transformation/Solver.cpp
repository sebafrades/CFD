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

constexpr int N = 64;
constexpr int M = 16;

void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename);

int main(){

    Matrix x = createMatrix(M+1, N+1, 'x');
    Matrix y = createMatrix(M+1, N+1, 'y');

    for (int j = 0; j < N+1; j++) {
        cout << x[M][j] << " ";
    }
    cout << endl;

    iterateMatrix(x, y, 1e-5);

    saveMatrix(x, y, "grid.csv");
}

// function definition
Matrix createMatrix(int rows, int cols, char type) {
    
    Matrix matrix(rows, vector<double>(cols, 0.0));

    switch(type){
        case 'x': 

            // superior and inferior row boundary
            for(int i = 0; i < cols; i++){

                if (i <= 16){  // first 17 points (0–16)

                    int idx = 17 - i;  // maps: i=0→17, i=16→1
                    double s = algebraicFunction(idx);

                    matrix[0][i] = 1.0 * (1.0 - s);
                    matrix[rows-1][i] = matrix[0][i];
                }

                else if (i <= 48){  // middle region (17–48)

                    matrix[0][i] = 1.0 + (1.0)*(double(i - 16)/32.0);
                    matrix[rows-1][i] = matrix[0][i];
                }

                else{  // last 17 points (49–64)

                    int idx = i - 47;  // maps: i=48→1, i=64→17
                    double s = algebraicFunction(idx);

                    matrix[0][i] = (1.0 + 0.5 + 0.5) + 1.0 * s;
                    matrix[rows-1][i] = matrix[0][i];
                }
            }

            // left and right column boundary
            for(int i=0; i<rows; i++){
                matrix[i][0] = 0.0;
                matrix[i][cols-1] = (1.0+0.5+0.5) + 1.0*algebraicFunction(cols-48);
            }

            break;

        case 'y':
            // y boundary
            for(int i=0; i<cols; i++){
                matrix[0][i] = 1.0;

                if (i <= 16) {
                    matrix[rows-1][i] = 0.0;
                }
                else if (i <= 48){
                    double xi = (i - 17.0) / (48.0 - 17.0); // [0,1]
                    xi = 2.0*xi - 1.0;                      // [-1,1]

                    double H = 0.2;

                    matrix[rows-1][i] = H * (1.0 - xi*xi);
                }
                else{
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

void iterateMatrix(Matrix& x, Matrix& y,double error) {

    double residual = 1.0;

    int rows = x.size();
    int cols = x[0].size();

    while (residual > error){

        residual = 0.0;

        for (int i = 1; i < rows-1; i++){
            for (int j = 1; j < cols-1; j++){
                
                double x_eta = (x[i+1][j] - x[i-1][j])/2;
                double x_xi = (x[i][j+1] - x[i][j-1])/2;
                double y_eta = (y[i+1][j] - y[i-1][j])/2;
                double y_xi = (y[i][j+1] - y[i][j-1])/2;

                double x_xi_eta = 0.25*(x[i+1][j+1] - x[i+1][j-1] - x[i-1][j+1] + x[i-1][j-1]);
                double y_xi_eta = 0.25*(y[i+1][j+1] - y[i+1][j-1] - y[i-1][j+1] + y[i-1][j-1]);

                double alpha = (x_eta)*(x_eta) + (y_eta)*(y_eta);
                double beta = (x_xi*x_eta) + (y_xi*y_eta);
                double gamma = (x_xi)*(x_xi) + (y_xi)*(y_xi);

                double x_new_val = (2*beta*x_xi_eta - alpha*(x[i][j+1]+x[i][j-1]) - gamma*(x[i+1][j]+x[i-1][j]))/(-2*(alpha+gamma));
                double y_new_val = (2*beta*y_xi_eta - alpha*(y[i][j+1]+y[i][j-1]) - gamma*(y[i+1][j]+y[i-1][j]))/(-2*(alpha+gamma));

                double diff = abs(x_new_val - x[i][j]);

                if (diff > residual){
                    residual = diff;
                }

                x[i][j] = x_new_val;
                y[i][j] = y_new_val;
            }
        }

        //cout << "Residual: " << residual << endl;
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


