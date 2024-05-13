#include "fdtd3d.h"
#include <iostream>
#include <fstream>

void update_Er(double ***Er, double ****Hth, double ****Hph, double ****check, int n){
    int OLD = (n + 1) % 2;
    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L -1 ; k++){
                Er[i][j][k] = Er[i][j][k] + dt / EPS0 / dth / r(i + 0.5) / std::sin(theta(j))
                                * ( std::sin(theta(j + 0.5)) * Hph[OLD][i][j][k] - std::sin(theta(j - 0.5)) * Hph[OLD][i][j-1][k] )
                                - dt / EPS0 / dph / r(i + 0.5) / std::sin(theta(j)) * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
            }
        }
    }
}

void update_Eth(double ****Eth, double ***Hr, double ****Hph, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = PML_L; j < Nth - PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + dt / EPS0 / dph / r(i) / std::sin(theta(j+0.5)) * ( Hr[i][j][k] - Hr[i][j][k-1])
                               - dt / EPS0 / dr / r(i) * ( r(i + 0.5) * Hph[OLD][i][j][k] - r(i - 0.5) * Hph[OLD][i-1][j][k] );
            }
        }
    }
}

void update_Eph(double ****Eph, double ***Hr, double ****Hth, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + dt / EPS0 / dr / r(i) * (r(i + 0.5) * Hth[OLD][i][j][k] - r(i - 0.5) * Hth[OLD][i-1][j][k])
                                - dt / EPS0 / dth / r(i) * (Hr[i][j][k] - Hr[i][j-1][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}