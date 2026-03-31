#include <iostream>

using namespace std;

#define N 16
#define M 8

int main(){
    double x[M+1][N+1];
    double y[M+1][N+1];

    // valores de x e y iniciales
    for(int i=0; i<M+1; i++){
        for(int j=0; j<N+1; j++){
            x[i][j] = 0.0;
            y[i][j] = 0.0;
        }
    }

    // contorno de x
    // fila superior e inferior
    for(int i=0; i<N+1; i++){
        x[0][i] = 0.0 + double(i);
        x[M][i] = 0.0 + double(i);
    }

    // columna izquierda y derecha
    for(int i=0; i<M+1; i++){
        x[i][0] = 0.0;
        x[i][N] = double(N);
    }

    // contorno de y
    // fila superior e inferior
    for(int i=0; i<N+1; i++){
        y[0][i] = double(M);

        if (i < N/4) {
            y[M][i] = 0.0;
        }
        else if (i >= N/4 && i < 3*N/8){
            y[M][i] = double(i - N/4);
        }
        else if (i >= 3*N/8 && i < 5*N/8){
            y[M][i] = double(2*M/8);
        }
        else if (i >= 5*N/8 && i < 3*N/4){
            y[M][i] = double(2*M/8 - (i - 5*N/8));
        }
        else if (i >= 3*N/4 && i < N){
            y[M][i] = 0.0;
        }
    }

    // columna izquierda y derecha
    for(int i=0; i<M+1; i++){
        y[i][0] = double(M-i);
        y[i][N] = double(M-i);
    } 
    
    for(int i=0; i<M+1; i++){
        for(int j=0; j<N+1; j++){
            cout << x[i][j] << " ";
        }
        cout << endl;
    }
}



