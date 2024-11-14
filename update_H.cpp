#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <omp.h>

void update_Hr(double ***Hr, double ****Eth, double ****Eph, double ****check, int n){
    int NEW = n % 2;
    omp_set_num_threads(10);
    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = PML_L; j < Nth - PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hr[i][j][k] = Hr[i][j][k] - dt / MU0 / dth / r(i) / std::sin(theta(j+0.5)) 
                                * (std::sin(theta(j+1.0)) * Eph[NEW][i][j+1][k] - std::sin(theta(j)) * Eph[NEW][i][j][k])
                                + dt / MU0 / dph / r(i) / std::sin(theta(j + 0.5)) * ( Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
                // Hr[PML_L+1][j][k] = 0.0;
            }
        }
    }
}

void update_Hth(double ****Hth, double ***Er, double ****Eph, double **Rs, double **Ls, double ***Bth, double **Bthr, double **Bthph, double *CHTHPH_00, double *CHTHPH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    // double c1 = r(1.0) * dt / r(0.5) / dr;
    #pragma omp parallel for collapse(2)
    for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
        for(int k = PML_L; k < Nph - PML_L; k++){
            // double c1 = ( MU0 * r(0.5) * dr - dt * r(0.0) * (Rs[j][k] / 2.0 - Ls[j][k] / dt) )
            //             / ( MU0 * r(0.5) * dr + dt * r(0.0) * (Rs[j][k] / 2.0 + Ls[j][k] / dt) );

            // double c2 = dt * r(1.0) / ( MU0 * r(0.5) * dr + dt * r(0.0) * (Rs[j][k] / 2.0 + Ls[j][k] / dt) );

            // Hth[NEW][0][j][k] = c1 * Hth[OLD][0][j][k] + c2 * Eph[NEW][1][j][k];


            double c1 = ( MU0 * r(0.5) * dr - dt * r(0.0) * (Rs[j][k] / 2.0 - Ls[j][k] / dt) )
                        / ( MU0 * r(0.5) * dr + dt * r(0.0) * (Rs[j][k] / 2.0 + Ls[j][k] / dt) );

            double c2 = dt / (MU0 * r(0.5) + dt / dr * r(0.0) * (Rs[j][k] / 2.0 + Ls[j][k] / dt));
            Hth[NEW][0][j][k] = c1 * Hth[OLD][0][j][k] 
                            - c2 * ( (Er[0][j][k+1] - Er[0][j][k]) / dph / std::sin(theta(j)) - Eph[NEW][1][j][k] * r(1.0) / dr );
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr - PML_L; i++){
        for(int j = PML_L + 1; j <= Nth - PML_L - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth[NEW][i][j][k] = Hth[OLD][i][j][k] - dt / MU0 / dph / r(i + 0.5) / std::sin(theta(j)) * ( Er[i][j][k+1] - Er[i][j][k] )
                                + dt / MU0 / dr / r(i + 0.5) * ( r(i + 1.0) * Eph[NEW][i+1][j][k] - r(i) * Eph[NEW][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Hph(double ****Hph, double ***Er, double ****Eth, double **Rs, double **Ls, double ***Bph, double **Bphr, double **Bphth, double *CHPHTH_00, double *CHPHTH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    // double c1 = r(1.0) * dt / r(0.5) / dr;
    #pragma omp parallel for collapse(2)
    for(int j = PML_L; j < Nth - PML_L; j++){
        for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
            double c1 = (MU0 * r(0.5) * dr - dt * r(0.0) * (Rs[j][k] / 2.0 - Ls[j][k] / dt)) 
                        / (MU0 * r(0.5) * dr + dt * r(0.0) * (Rs[j][k] / 2.0 + Ls[j][k] / dt));
            
            double c2 = dt / (MU0 * r(0.5) + dt / dr * r(0.0) * (Rs[j][k] / 2.0 + Ls[j][k] / dt));
            
            Hph[NEW][0][j][k] = c1 * Hph[OLD][0][j][k] + 
                                        c2 * ( - r(1.0) / dr * Eth[NEW][1][j][k] + ( Er[0][j+1][k] - Er[0][j][k] ) / dth);
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr - PML_L; i++){
        for(int j = PML_L; j < Nth - PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph[NEW][i][j][k] = Hph[OLD][i][j][k] - dt / MU0 / dr / r(i + 0.5) * ( r(i + 1.0) * Eth[NEW][i+1][j][k] - r(i) * Eth[NEW][i][j][k] )
                                + dt / MU0 / dth / r(i+0.5) * ( Er[i][j+1][k] - Er[i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}