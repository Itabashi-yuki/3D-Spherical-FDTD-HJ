#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <string>

void update_Dr_PML(double ****Drth1, double ****Drth2, double ****Drph, double ****Dr, double ***Hr, double ****Hth, double ****Hph,
                     double *CDRTH1_00, double *CDRTH1_01, double *CDRPH_00, double *CDRPH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    // std::ofstream ofs_Dr_check(PATH + "data/" + global_dirName + "Dr_check.dat",std::ios::app);

    for(int i = 0; i < Nr; i++){    /* Dr(i+1/2, j, k) */
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 1; k <= Nph - 1; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];

                // if(i == 0 && j == 50 && k < 13){
                //     ofs_Dr_check << n << " " << k * R0 * dph * 1e-3 << " " 
                //                     << Drth1[NEW][i][j][k] << " " << Drth2[NEW][i][j][k] << " " << Drph[NEW][i][j][k] << " " 
                //                     << Hph[NEW][i][j][k] << " " << Hph[NEW][i][j-1][k] << " " << Hph[OLD][i][j][k] << " "  << Hth[OLD][i][j][k] << " " 
                //                     << Hth[OLD][i][j-1][k] << std::endl;
                // }    
                // check[NEW][i][j][k] += 1.0;
                                //     ofs_Dr_check << n << " " << i << " " << j << " " << k << " " << Drth1[NEW][i][j][k] << " " << Drth2[NEW][i][j][k] << " " << Drph[NEW][i][j][k] << " " << Dr[NEW][i][j][k] << " "
                                // << Hph[OLD][i][j][k] << " " << Hph[OLD][i][j-1][k] << " " << Hth[OLD][i][j][k] << " " << Hth[OLD][i][j][k-1] <<  std::endl;
                
            }
        }
    }

}

void update_Dth_PML(double ****Dthph, double ****Dthr, double ****Dthr_tilde, double ****Dth, double ***Hr, double ***Hph_tilde, double *CDTHPH_00, double *CDTHPH_01,
                     double *CDTHR_10, double *CDTHR_11, double *CDTHR_TILDE_00, double *CDTHR_TILDE_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n+1) % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= Nph - 1; k++){
                Dthph[NEW][i][j][k] = CDTHPH_00[k] * Dthph[OLD][i][j][k] + CDTHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph
                                     * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Dthr_tilde[NEW][i][j][k] = CDTHR_TILDE_00[i] * Dthr_tilde[OLD][i][j][k] - CDTHR_TILDE_01[i] / dr
                                     * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Dthr[NEW][i][j][k] = CDTHR_10[i] * Dthr[OLD][i][j][k] + CDTHR_11[i] / dt
                                     * ( Dthr_tilde[NEW][i][j][k] - Dthr_tilde[OLD][i][j][k] );
                Dth[NEW][i][j][k] = Dthph[NEW][i][j][k] + Dthr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;

            }
        }
    }
}

void update_Dph_PML(double ****Dphr, double ****Dphr_tilde, double ****Dphth, double ****Dph, double ***Hr, double ***Hth_tilde,
                     double *CDPHR_10, double *CDPHR_11, double *CDPHR_TILDE_00, double *CDPHR_TILDE_01, double *CDPHTH_00,
                     double *CDPHTH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < Nph; k++){
                Dphth[NEW][i][j][k] = CDPHTH_00[j] * Dphth[OLD][i][j][k] - CDPHTH_01[j] / r(i) / dth * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Dphr_tilde[NEW][i][j][k] = CDPHR_TILDE_00[i] * Dphr_tilde[OLD][i][j][k] + CDPHR_TILDE_01[i] / dr
                                         * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Dphr[NEW][i][j][k] = CDPHR_10[i] * Dphr[OLD][i][j][k] + CDPHR_11[i] / dt * ( Dphr_tilde[NEW][i][j][k] - Dphr_tilde[OLD][i][j][k] );
                Dph[NEW][i][j][k] = Dphr[NEW][i][j][k] + Dphth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}