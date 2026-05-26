#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

using Matrix = vector<vector<double>>;

// matrix creation function declaration
Matrix createMatrix(int rows, int cols, char type);

// matrix iteration process function declaration
void iterateMatrix(Matrix& x, Matrix& y,double error);

constexpr int N = 32;
constexpr int M = 16;

void saveMatrix(const Matrix& x, const Matrix& y, const std::string& filename);

int main(){

    Matrix x = createMatrix(M+1, N+1, 'x');
    Matrix y = createMatrix(M+1, N+1, 'y');

    iterateMatrix(x, y, 1e-5);

    saveMatrix(x, y, "grid.csv");
}

// function definition
Matrix createMatrix(int rows, int cols, char type) {
    
    Matrix matrix(rows, vector<double>(cols, 0.0));

    // valores iniciales
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            matrix[i][j] = 0.0;
        }
    }

    switch(type){
        case 'x': 

            // fila superior e inferior
            for(int i=0; i<cols; i++){
                matrix[0][i] = 0.0 + double(i);
                matrix[rows-1][i] = 0.0 + double(i);
            }

            // columna izquierda y derecha
            for(int i=0; i<rows; i++){
                matrix[i][0] = 0.0;
                matrix[i][cols-1] = double(cols-1);
            }

            break;

        case 'y':
            // contorno de y
            // fila superior e inferior
            for(int i=0; i<cols; i++){
                matrix[0][i] = double(M);

                if (i < N/4) {
                    matrix[rows-1][i] = 0.0;
                }
                else if (i >= N/4 && i < 3*N/8){
                    matrix[rows-1][i] = double(i - N/4);
                }
                else if (i >= 3*N/8 && i < 5*N/8){
                    matrix[rows-1][i] = double(2*M/8);
                }
                else if (i >= 5*N/8 && i < 3*N/4){
                    matrix[rows-1][i] = double(2*M/8 - (i - 5*N/8));
                }
                else if (i >= 3*N/4 && i < N){
                    matrix[rows-1][i] = 0.0;
                }
            }

            // columna izquierda y derecha
            for(int i=0; i<rows; i++){
                matrix[i][0] = double(M-i);
                matrix[i][cols-1] = double(M-i);
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

