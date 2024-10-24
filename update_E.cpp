#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <omp.h>

void update_Er(double ***Er, double ****Hth, double ****Hph, double ****Jr, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;    
    omp_set_num_threads(8);
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L; i++){
        for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L -1 ; k++){
                Er[i][j][k] = Er[i][j][k] + dt / EPS0 / dth / r(i + 0.5) / std::sin(theta(j))
                                * ( std::sin(theta(j + 0.5)) * Hph[OLD][i][j][k] - std::sin(theta(j - 0.5)) * Hph[OLD][i][j-1][k] )
                                - dt / EPS0 / dph / r(i + 0.5) / std::sin(theta(j)) * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] ) - dt / EPS0 * Jr[OLD][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
    
    // #pragma omp parallel for collapse(3)
    // for(int i = 0; i < PML_L; i++){
    //     for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
    //         for(int k = PML_L + 1; k <= Nph - PML_L -1 ; k++){
    //             Er[i][j][k] = 0.0;
    //             // check[NEW][i][j][k] += 2.0;
    //         }
    //     }
    // }
}

void update_Eth(double ****Eth, double ***Hr, double ****Hph, double ****Jth, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    // omp_set_num_threads(10);
    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = PML_L; j < Nth - PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + dt / EPS0 / dph / r(i) / std::sin(theta(j+0.5)) * ( Hr[i][j][k] - Hr[i][j][k-1])
                               - dt / EPS0 / dr / r(i) * ( r(i + 0.5) * Hph[OLD][i][j][k] - r(i - 0.5) * Hph[OLD][i-1][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Eth[NEW][PML_L + 1][j][k] = 0.0;
            }
        }
    }
    
    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i <= PML_L; i++){
    //     for(int j = PML_L; j < Nth - PML_L; j++){
    //         for(int k = PML_L + 1; k <= Nph - PML_L -1 ; k++){
    //             Eth[NEW][i][j][k] = 0.0;
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }
}

void update_Eph(double ****Eph, double ***Hr, double ****Hth, double ****Jph, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;
    
    // omp_set_num_threads(10);
    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + dt / EPS0 / dr / r(i) * (r(i + 0.5) * Hth[OLD][i][j][k] - r(i - 0.5) * Hth[OLD][i-1][j][k])
                                - dt / EPS0 / dth / r(i) * (Hr[i][j][k] - Hr[i][j-1][k]) - dt / EPS0 * Jph[OLD][i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Eph[NEW][PML_L + 1][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i <= PML_L; i++){
    //     for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
    //         for(int k = PML_L; k < Nph - PML_L ; k++){
    //             Eph[NEW][i][j][k] = 0.0;
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }
}