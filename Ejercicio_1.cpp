#include <iostream>

using namespace std;

#define N 16
#define M 8

int main(){
    double x[M][N];
    double y[M][N];

    // valores de x e y iniciales
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            x[i][j] = 0.0;
            y[i][j] = 0.0;
        }
    }

    // contorno de x
    // fila superior e inferior
    for(int i=0; i<N; i++){
        x[0][i] = 0.0 + double(i);
        x[M-1][i] = 0.0 + double(i);
    }

    // columna izquierda y derecha
    for(int i=0; i<M; i++){
        x[i][0] = 0.0;
        x[i][N-1] = double(N);
    }

    // contorno de y
    // fila superior e inferior
    for(int i=0; i<N; i++){
        y[0][i] = double(M);

        if (i < N/4) {
            y[M-1][i] = 0.0
        }
        else if (i >= N/4 && i < 3*N/8){

        }
        else if (i >= 3*N/8 && i < N/2){

        }
        else if (i >= N/2 && i < 5*N/8){

        }
        else if (i >= 5*N/8 && i < 3*N/4){

        }
        else if (i >= 3*N/4 && i < 7*N/8){

        }
        else if (i >= 7*N/8 && i < N){

        }
        y[M-1][i] = ;
    }




    for(int i=0; i<N; i++){
        cout << x[0][i] << " ";
    }        
}



