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

    for(int i=0; i<N; i++){
        x[0][i] = 0.0 + double(i);
        x[M-1][i] = 0.0 + double(i);
    }

    for(int i=0; i<N; i++){
        cout << x[0][i] << " ";
    }        
}



