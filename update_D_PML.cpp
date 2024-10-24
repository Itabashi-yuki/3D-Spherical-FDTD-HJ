#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <string>

void update_Dr_PML(double ****Drth1, double ****Drth2, double ****Drph, double ****Dr, double ***Hr, double ****Hth, double ****Hph,
                     double *CDRTH1_00, double *CDRTH1_01, double *CDRPH_00, double *CDRPH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    omp_set_num_threads(8);
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 1; k <= PML_L; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;      
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;   
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;   
            }
        }
    }
                // if(n >= 39 && n <= 42){
                //     std::ofstream ofs_PML("./data/" + global_dirName + "/PML_" + std::to_string(n) +".dat",std::ios::app);
                //     for(int i = 84; i < Nr; i++){
                //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                //             ofs_PML << i * dr * 1.0e-3 << " " << k * R0 * dph * 1.0e-3 << " " << Erth1[NEW][i][Nth/2][k] << " " << Erth2[NEW][i][Nth/2][k] 
                //                     << " " << Erph[NEW][i][Nth/2][k] << " " << Er[i][Nth/2][k] << " " << Hph[OLD][i][Nth/2][k] << " " << Hph[OLD][i][Nth/2 - 1][k] << std::endl;
                //         }
                //         ofs_PML << std::endl;
                //     }
                // }

    #pragma omp parallel for collapse(3)
    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;   
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Drth1[NEW][i][j][k] = CDRTH1_00[j] * Drth1[OLD][i][j][k] + CDRTH1_01[j] / r(i + 0.5) / dth  
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Drth2[NEW][i][j][k] = Drth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Drph[NEW][i][j][k] = CDRPH_00[k] * Drph[OLD][i][j][k] - CDRPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Dr[NEW][i][j][k] = Drth1[NEW][i][j][k] + Drth2[NEW][i][j][k] + Drph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;   
            }
        }
    }   
    // exit(0);
}

void update_Dth_PML(double ****Dthph, double ****Dthr, double ****Dthr_tilde, double ****Dth, double ***Hr, double ***Hph_tilde, double *CDTHPH_00, double *CDTHPH_01,
                     double *CDTHR_10, double *CDTHR_11, double *CDTHR_TILDE_00, double *CDTHR_TILDE_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n+1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
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

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
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

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= PML_L; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
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

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
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

    #pragma omp parallel for collapse(3)
    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
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

    #pragma omp parallel for collapse(3)
    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
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

                //     if(n >= 39 && n <= 42){
                //     std::ofstream ofs_PML("./data/" + global_dirName + "/PML_" + std::to_string(n) +".dat",std::ios::app);
                //     for(int i = 0; i < Nr; i++){
                //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                //             ofs_PML << i * dr * 1.0e-3 << " " << k * R0 * dph * 1.0e-3 << " " <<  Ethr_tilde[NEW][i][Nth/2][k] << std::endl;
                //         }
                //         ofs_PML << std::endl;
                //     }
                // }
}

void update_Dph_PML(double ****Dphr, double ****Dphr_tilde, double ****Dphth, double ****Dph, double ***Hr, double ***Hth_tilde,
                     double *CDPHR_10, double *CDPHR_11, double *CDPHR_TILDE_00, double *CDPHR_TILDE_01, double *CDPHTH_00,
                     double *CDPHTH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Dphth[NEW][i][j][k] = CDPHTH_00[j] * Dphth[OLD][i][j][k] - CDPHTH_01[j] / r(i) / dth * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Dphr_tilde[NEW][i][j][k] = CDPHR_TILDE_00[i] * Dphr_tilde[OLD][i][j][k] + CDPHR_TILDE_01[i] / dr
                                         * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Dphr[NEW][i][j][k] = CDPHR_10[i] * Dphr[OLD][i][j][k] + CDPHR_11[i] / dt * ( Dphr_tilde[NEW][i][j][k] - Dphr_tilde[OLD][i][j][k] );
                Dph[NEW][i][j][k] = Dphr[NEW][i][j][k] + Dphth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Dphth[NEW][i][j][k] = CDPHTH_00[j] * Dphth[OLD][i][j][k] - CDPHTH_01[j] / r(i) / dth * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Dphr_tilde[NEW][i][j][k] = CDPHR_TILDE_00[i] * Dphr_tilde[OLD][i][j][k] + CDPHR_TILDE_01[i] / dr
                                         * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Dphr[NEW][i][j][k] = CDPHR_10[i] * Dphr[OLD][i][j][k] + CDPHR_11[i] / dt * ( Dphr_tilde[NEW][i][j][k] - Dphr_tilde[OLD][i][j][k] );
                Dph[NEW][i][j][k] = Dphr[NEW][i][j][k] + Dphth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Dphth[NEW][i][j][k] = CDPHTH_00[j] * Dphth[OLD][i][j][k] - CDPHTH_01[j] / r(i) / dth * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Dphr_tilde[NEW][i][j][k] = CDPHR_TILDE_00[i] * Dphr_tilde[OLD][i][j][k] + CDPHR_TILDE_01[i] / dr
                                         * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Dphr[NEW][i][j][k] = CDPHR_10[i] * Dphr[OLD][i][j][k] + CDPHR_11[i] / dt * ( Dphr_tilde[NEW][i][j][k] - Dphr_tilde[OLD][i][j][k] );
                Dph[NEW][i][j][k] = Dphr[NEW][i][j][k] + Dphth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Dphth[NEW][i][j][k] = CDPHTH_00[j] * Dphth[OLD][i][j][k] - CDPHTH_01[j] / r(i) / dth * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Dphr_tilde[NEW][i][j][k] = CDPHR_TILDE_00[i] * Dphr_tilde[OLD][i][j][k] + CDPHR_TILDE_01[i] / dr
                                         * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Dphr[NEW][i][j][k] = CDPHR_10[i] * Dphr[OLD][i][j][k] + CDPHR_11[i] / dt * ( Dphr_tilde[NEW][i][j][k] - Dphr_tilde[OLD][i][j][k] );
                Dph[NEW][i][j][k] = Dphr[NEW][i][j][k] + Dphth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Dphth[NEW][i][j][k] = CDPHTH_00[j] * Dphth[OLD][i][j][k] - CDPHTH_01[j] / r(i) / dth * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Dphr_tilde[NEW][i][j][k] = CDPHR_TILDE_00[i] * Dphr_tilde[OLD][i][j][k] + CDPHR_TILDE_01[i] / dr
                                         * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Dphr[NEW][i][j][k] = CDPHR_10[i] * Dphr[OLD][i][j][k] + CDPHR_11[i] / dt * ( Dphr_tilde[NEW][i][j][k] - Dphr_tilde[OLD][i][j][k] );
                Dph[NEW][i][j][k] = Dphr[NEW][i][j][k] + Dphth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
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