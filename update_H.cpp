#include "fdtd3d.h"
#include <iostream>
#include <fstream>

void update_Hr(double ***Hr, double ****Eth, double ****Eph, double ****check, int n){
    int NEW = n % 2;

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = PML_L; j < Nth - PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hr[i][j][k] = Hr[i][j][k] - dt / MU0 / dth / r(i) / std::sin(theta(j+0.5)) 
                                * (std::sin(theta(j+1.0)) * Eph[NEW][i][j+1][k] - std::sin(theta(j)) * Eph[NEW][i][j][k])
                                + dt / MU0 / dph / r(i) / std::sin(theta(j + 0.5)) * ( Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Hth(double ****Hth, double ***Er, double ****Eph, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth[NEW][i][j][k] = Hth[OLD][i][j][k] - dt / MU0 / dph / r(i + 0.5) / std::sin(theta(j)) * ( Er[i][j][k+1] - Er[i][j][k] )
                                + dt / MU0 / dr / r(i + 0.5) * ( r(i + 1.0) * Eph[NEW][i+1][j][k] - r(i) * Eph[NEW][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Hph(double ****Hph, double ***Er, double ****Eth, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = PML_L; j < Nth - PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph[NEW][i][j][k] = Hph[OLD][i][j][k] - dt / MU0 / dr / r(i + 0.5) * ( r(i + 1.0) * Eth[NEW][i+1][j][k] - r(i) * Eth[NEW][i][j][k] )
                                + dt / MU0 / dth / r(i+0.5) * ( Er[i][j+1][k] - Er[i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}