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

constexpr int N = 30, M = 30;

const int H = 1;
double a_ang = 26.565; // angle in degrees

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
        case 'x': {

            double tan_a = tan(a_ang * M_PI / 180.0);

            // top & bottom (ξ direction)
            for(int j = 0; j < cols; j++){
                double s = double(j) / (cols - 1);

                matrix[0][j] = s * H * (1 + tan_a);              // top
                matrix[rows-1][j] = H * tan_a + s * H;           // bottom
            }

            // left & right (η direction)
            for(int i = 0; i < rows; i++){
                double t = double(i) / (rows - 1);

                matrix[i][0] = t * H * tan_a;                    // left
                matrix[i][cols-1] = H * (1 + tan_a);             // right (vertical line)
            }

            break;
        }

        case 'y': {

            // top & bottom
            for(int j = 0; j < cols; j++){
                matrix[0][j] = H;
                matrix[rows-1][j] = 0.0;
            }

            // left & right
            for(int i = 0; i < rows; i++){
                double t = double(i) / (rows - 1);

                matrix[i][0] = H * (1.0 - t);
                matrix[i][cols-1] = H * (1.0 - t);
            }

            break;
        }
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

                double J = x_xi*y_eta - x_eta*y_xi;
                double P = (2*x_eta*x_xi_eta)/(x_eta*x_eta + y_eta*y_eta);

                double x_new_val = (-J*J*P*x_xi + 2*beta*x_xi_eta - alpha*(x[i][j+1]+x[i][j-1]) - gamma*(x[i+1][j]+x[i-1][j]))/(-2*(alpha+gamma));
                double y_new_val = (-J*J*P*x_xi + 2*beta*y_xi_eta - alpha*(y[i][j+1]+y[i][j-1]) - gamma*(y[i+1][j]+y[i-1][j]))/(-2*(alpha+gamma));

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


