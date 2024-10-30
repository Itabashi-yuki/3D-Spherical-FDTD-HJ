#include "fdtd3d.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

void update_Hr_PML(double ***Hr, double ***Hrth1, double ***Hrth2, double ***Hrph,
                     double ****Eth, double ****Eph, double *CHRTH1_00, double *CHRTH1_01,
                      double *CHRPH_00, double *CHRPH_01, double ****check, int n){
    int NEW = n % 2;
    
    omp_set_num_threads(8);
    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 0; k < PML_L; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Hr[PML_L+1][j][k] = 0.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Hr[PML_L+1][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i <= PML_L; i++){
    //     for(int j = 0; j < Nth; j++){
    //         for(int k = PML_L; k < Nph - PML_L; k++){
    //             Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
    //             Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
    //             Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
    //             Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
           }
        }
    }
}


void update_Hth_PML(double ****Hth, double ***Hthr, double ***Hthph, double ****Hthr_tilde, double ***Er, double ****Eph, double ****Eph_tilde,
                 double ***Bth, double **Bthr, double **Bthph, double **Rs, double **Ls,
                 double *CHTHPH_00,double *CHTHPH_01, double *CHTHR_TILDE_00, double *CHTHR_TILDE_01,
                double *CHTHR_10, double *CHTHR_11, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    double c1 = r(1.0) * dt / r(0.5) / dr;
    for(int j = 1; j <= Nth - 1; j++){
        for(int k = 0; k < PML_L; k++){
            Bthr[j][k] = Bthr[j][k] + c1 * Eph[NEW][1][j][k];
            Bthph[j][k] = CHTHPH_00[k] * Bthph[j][k] - CHTHPH_01[k] / r(0.5) / std::sin(theta(j)) / dph * 
                            (Er[0][j][k+1] - Er[0][j][k]);
            Bth[NEW][j][k] = Bthr[j][k] + Bthph[j][k];

            // double alpha = r(0.0) * Rs[j][k] - MU0;
            // double beta = r(0.0) * Rs[j][k];
            // double ch1 = (alpha / dt - beta / 2.0) / (alpha / dt + beta / 2.0);
            // double ch2 = 1.0 / (alpha / dt + beta / 2.0) / dt;

            double alpha = (r(0.0) * Rs[j][k]) / (r(0.5) * dr);
            double beta = (r(0.0) * Ls[j][k]) / (r(0.5) * dr) + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            // std::cout << r(0.0) << " " << Ls[j][k] << " " << Rs[j][k] << " " <<  MU0 << std::endl;
            // std::cout << alpha << " " << beta << " " << ch1 << " " << ch2 <<  std::endl;
            // exit(0);
            Hth[NEW][0][j][k] = ch1 * Hth[OLD][0][j][k] + ch2 * (Bth[NEW][j][k] - Bth[OLD][j][k]); 
            // check[NEW][0][j][k] += 1.0;
        }
    }

    for(int j = 1; j <= Nth - 1; j++){
        for(int k = Nph - PML_L; k < Nph; k++){
            Bthr[j][k] = Bthr[j][k] + c1 * Eph[NEW][1][j][k];
            Bthph[j][k] = CHTHPH_00[k] * Bthph[j][k] - CHTHPH_01[k] / r(0.5) / std::sin(theta(j)) / dph * 
                            (Er[0][j][k+1] - Er[0][j][k]);
            Bth[NEW][j][k] = Bthr[j][k] + Bthph[j][k];

            double alpha = (r(0.0) * Rs[j][k]) / (r(0.5) * dr);
            double beta = (r(0.0) * Ls[j][k]) / (r(0.5) * dr) + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            Hth[NEW][0][j][k] = ch1 * Hth[OLD][0][j][k] + ch2 * (Bth[NEW][j][k] - Bth[OLD][j][k]); 
            // check[NEW][0][j][k] += 1.0;
        }
    }

    for(int j = 1; j <= PML_L; j++){
        for(int k = PML_L; k < Nph - PML_L; k++){
            Bthr[j][k] = Bthr[j][k] + c1 * Eph[NEW][1][j][k];
            Bthph[j][k] = CHTHPH_00[k] * Bthph[j][k] - CHTHPH_01[k] / r(0.5) / std::sin(theta(j)) / dph * 
                            (Er[0][j][k+1] - Er[0][j][k]);
            Bth[NEW][j][k] = Bthr[j][k] + Bthph[j][k];

            double alpha = (r(0.0) * Rs[j][k]) / (r(0.5) * dr);
            double beta = (r(0.0) * Ls[j][k]) / (r(0.5) * dr) + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            Hth[NEW][0][j][k] = ch1 * Hth[OLD][0][j][k] + ch2 * (Bth[NEW][j][k] - Bth[OLD][j][k]); 
            // check[NEW][0][j][k] += 1.0;
        }
    }

    for(int j = Nth - PML_L; j <= Nth - 1; j++){
        for(int k = PML_L; k < Nph - PML_L; k++){
            Bthr[j][k] = Bthr[j][k] + c1 * Eph[NEW][1][j][k];
            Bthph[j][k] = CHTHPH_00[k] * Bthph[j][k] - CHTHPH_01[k] / r(0.5) / std::sin(theta(j)) / dph * 
                            (Er[0][j][k+1] - Er[0][j][k]);
            Bth[NEW][j][k] = Bthr[j][k] + Bthph[j][k];

            double alpha = (r(0.0) * Rs[j][k]) / (r(0.5) * dr);
            double beta = (r(0.0) * Ls[j][k]) / (r(0.5) * dr) + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            Hth[NEW][0][j][k] = ch1 * Hth[OLD][0][j][k] + ch2 * (Bth[NEW][j][k] - Bth[OLD][j][k]); 
            // check[NEW][0][j][k] += 1.0;
        }
    }

    std::ofstream ofs_Hth(PATH + "data/" + global_dirName + "Hth_check.dat",std::ios::app);

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] 
                                - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] 
                                        + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] 
                                + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];

                // if(i == 1 && j == 50){
                //     ofs_Hth << n << " " << k << " " << Hthph[i][j][k] << " " << Hthr_tilde[NEW][i][j][k] << " " << Hthr[i][j][k] << " " << Hth[NEW][i][j][k] << " "
                //             << CHTHPH_00[k] << " " << CHTHPH_01[k] << " " << CHTHR_TILDE_00[i] << " " << CHTHR_TILDE_01[i] << " " << CHTHR_10[i] << " " << CHTHR_11[i] << std::endl; 
                // }
                // check[NEW][i][j][k] += 1.0;
                // Hth[NEW][PML_L][j][k] = 0.0;
            }
        }
    }
    // ofs_Hth << std::endl;

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Hth[NEW][PML_L][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i < PML_L; i++){
    //     for(int j = 1; j <= Nth - 1; j++){
    //         for(int k = PML_L; k < Nph - PML_L; k++){
    //             Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
    //             Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
    //             Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
    //             Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    // #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];

                if(j == 50 && k == 50){
                    ofs_Hth << n << " " << i << " " << Hthph[i][j][k] << " " << Hthr_tilde[NEW][i][j][k] << " " << Hthr[i][j][k] << " " << Hth[NEW][i][j][k] << " "
                            << CHTHPH_00[k] << " " << CHTHPH_01[k] << " " << CHTHR_TILDE_00[i] << " " << CHTHR_TILDE_01[i] << " " << CHTHR_10[i] << " " << CHTHR_11[i] << std::endl; 
                }

                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr - PML_L; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;

                // if(i == 1 && k == 50){
                //     ofs_Hth << n << " " << j << " " << Hthph[i][j][k] << " " << Hthr_tilde[NEW][i][j][k] << " " << Hthr[i][j][k] << " " << Hth[NEW][i][j][k] << " "
                //             << CHTHPH_00[k] << " " << CHTHPH_01[k] << " " << CHTHR_TILDE_00[i] << " " << CHTHR_TILDE_01[i] << " " << CHTHR_10[i] << " " << CHTHR_11[i] << std::endl; 
                // }
            }
        }
    }
    ofs_Hth << std::endl;

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr - PML_L; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }   
}

void update_Hph_PML(double ****Hph, double ***Hphr, double ***Hphth, double ****Hphr_tilde, double ***Er, double ****Eth, double ****Eth_tilde,
                     double ***Bph, double **Bphr, double **Bphth, double **Rs, double **Ls,
                     double *CHPHTH_00, double *CHPHTH_01, double *CHPHR_TILDE_00, double *CHPHR_TILDE_01, double *CHPHR_10, double *CHPHR_11, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;
    // std::ofstream ofs_c(PATH + "data/" + global_dirName + "c1c2.dat",std::ios::app);

    double c1 = r(1.0) * dt / r(0.5) / dr;
    for(int j = 0; j < Nth; j++){
        for(int k = 1; k <= PML_L; k++){
            Bphr[j][k] = Bphr[j][k] - c1 * Eth[NEW][1][j][k];
            Bphth[j][k] = CHPHTH_00[j] * Bphth[j][k] + CHPHTH_01[j] / r(0.5) / dth * (Er[0][j+1][k] - Er[0][j][k]);
            Bph[NEW][j][k] = Bphr[j][k] + Bphth[j][k];

            double alpha = r(0.0) * Rs[j][k] / r(0.5) / dr;
            double beta = r(0.0) * Ls[j][k] / r(0.5) / dr + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            // double alpha = r(0.0) * (Rs[j][k] + Rs[j+1][k] + Rs[j][k+1] + Rs[j+1][k+1]) / 4.0 / r(0.5) / dr;
            // double beta = r(0.0) * (Ls[j][k] + Ls[j+1][k] + Ls[j][k+1] + Ls[j+1][k+1]) / 4.0 / r(0.5) / dr + MU0;
            // double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            // double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            // double alpha = r(0.0) * Rs[j][k] - MU0;
            // double beta = r(0.0) * Rs[j][k];
            // double ch1 = (alpha / dt - beta / 2.0) / (alpha / dt + beta / 2.0);
            // double ch2 = 1.0 / (alpha / dt + beta / 2.0) / dt;

            // std::cout << r(0.0) << " " << Ls[j][k] << " " << r(0.5) << " " << dr << " " << MU0 << std::endl;
            // std::cout << alpha << " " << beta << " " << ch1 << " " << ch2 <<  std::endl;
            // exit(0);
            Hph[NEW][0][j][k] = ch1 * Hph[OLD][0][j][k] + ch2 * (Bph[NEW][j][k] - Bph[OLD][j][k]);
            // check[NEW][0][j][k] += 1.0;
        }
    }
    // ofs_c << std::endl;

    for(int j = 0; j < Nth; j++){
        for(int k = Nph - PML_L; k <= Nph - 1; k++){
            Bphr[j][k] = Bphr[j][k] - c1 * Eth[NEW][1][j][k];
            Bphth[j][k] = CHPHTH_00[j] * Bphth[j][k] + CHPHTH_01[j] / r(0.5) / dth * (Er[0][j+1][k] - Er[0][j][k]);
            Bph[NEW][j][k] = Bphr[j][k] + Bphth[j][k];

            double alpha = r(0.0) * Rs[j][k] / r(0.5) / dr;
            double beta = r(0.0) * Ls[j][k] / r(0.5) / dr + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            // double alpha = r(0.0) * Rs[j][k] - MU0;
            // double beta = r(0.0) * Rs[j][k];
            // double ch1 = (alpha / dt - beta / 2.0) / (alpha / dt + beta / 2.0);
            // double ch2 = 1.0 / (alpha / dt + beta / 2.0) / dt;

            Hph[NEW][0][j][k] = ch1 * Hph[OLD][0][j][k] + ch2 * (Bph[NEW][j][k] - Bph[OLD][j][k]);
            // check[NEW][0][j][k] += 1.0;
        }
    }

    for(int j = 0; j < PML_L; j++){
        for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
            Bphr[j][k] = Bphr[j][k] - c1 * Eth[NEW][1][j][k];
            Bphth[j][k] = CHPHTH_00[j] * Bphth[j][k] + CHPHTH_01[j] / r(0.5) / dth * (Er[0][j+1][k] - Er[0][j][k]);
            Bph[NEW][j][k] = Bphr[j][k] + Bphth[j][k];

            double alpha = r(0.0) * Rs[j][k] / r(0.5) / dr;
            double beta = r(0.0) * Ls[j][k] / r(0.5) / dr + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            // double alpha = r(0.0) * Rs[j][k] - MU0;
            // double beta = r(0.0) * Rs[j][k];
            // double ch1 = (alpha / dt - beta / 2.0) / (alpha / dt + beta / 2.0);
            // double ch2 = 1.0 / (alpha / dt + beta / 2.0) / dt;

            Hph[NEW][0][j][k] = ch1 * Hph[OLD][0][j][k] + ch2 * (Bph[NEW][j][k] - Bph[OLD][j][k]);
            // check[NEW][0][j][k] += 1.0;
        }
    }

    for(int j = Nth - PML_L; j < Nth; j++){
        for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
            Bphr[j][k] = Bphr[j][k] - c1 * Eth[NEW][1][j][k];
            Bphth[j][k] = CHPHTH_00[j] * Bphth[j][k] + CHPHTH_01[j] / r(0.5) / dth * (Er[0][j+1][k] - Er[0][j][k]);
            Bph[NEW][j][k] = Bphr[j][k] + Bphth[j][k];

            double alpha = r(0.0) * Rs[j][k] / r(0.5) / dr;
            double beta = r(0.0) * Ls[j][k] / r(0.5) / dr + MU0;
            double ch1 = (beta / dt - alpha / 2.0) / (beta / dt + alpha / 2.0);
            double ch2 = 1.0 / (beta / dt + alpha / 2.0) / dt;

            // double alpha = r(0.0) * Rs[j][k] - MU0;
            // double beta = r(0.0) * Rs[j][k];
            // double ch1 = (alpha / dt - beta / 2.0) / (alpha / dt + beta / 2.0);
            // double ch2 = 1.0 / (alpha / dt + beta / 2.0) / dt;

            Hph[NEW][0][j][k] = ch1 * Hph[OLD][0][j][k] + ch2 * (Bph[NEW][j][k] - Bph[OLD][j][k]);
            // check[NEW][0][j][k] += 1.0;
        }
    }

    std::ofstream ofs_Hph(PATH + "data/" + global_dirName + "Hph_check.dat",std::ios::app);    

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Hph[NEW][PML_L][j][k] = 0.0;

                // if(i == 1 && j == 50){
                //     ofs_Hph << n << " " << k << " " << Hphth[i][j][k] << " " << Hphr_tilde[NEW][i][j][k] << " " << Hphr[i][j][k] << " " << Hph[NEW][i][j][k] << " "
                //             << CHPHTH_00[j] << " " << CHPHTH_01[j] << " " << CHPHR_TILDE_00[i] << " " << CHPHR_TILDE_01[i] << " " << CHPHR_10[i] << " " << CHPHR_11[i] << " " << std::endl;
                // } 
            }
        }
    }
    // ofs_Hph << std::endl;

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
                // Hph[NEW][PML_L][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i < PML_L; i++){
    //     for(int j = 0; j < Nth; j++){
    //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
    //             Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
    //             Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
    //             Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
    //             Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    // #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];

                if(j == 50 && k == 50){
                    ofs_Hph << n << " " << i << " " << Hphth[i][j][k] << " " << Hphr_tilde[NEW][i][j][k] << " " << Hphr[i][j][k] << " " << Hph[NEW][i][j][k] << " "
                            << CHPHTH_00[j] << " " << CHPHTH_01[j] << " " << CHPHR_TILDE_00[i] << " " << CHPHR_TILDE_01[i] << " " << CHPHR_10[i] << " " << CHPHR_11[i] << " " << std::endl;
                } 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

                // if(n >= 39 && n <= 42){
                //     std::ofstream ofs_PML("./data/" + global_dirName + "/PML_" + std::to_string(n) +".dat",std::ios::app);
                //     for(int i = 84; i < Nr; i++){
                //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                //             ofs_PML << i * dr * 1.0e-3 << " " << k * R0 * dph * 1.0e-3 << " " << Hphth[i][Nth/2][k] << " " << Hphr_tilde[NEW][i][Nth/2][k] 
                //             << " " << Hphr[i][Nth/2][k] << " " << Hph[NEW][i][Nth/2][k] << " " << Eth_tilde[NEW][i+1][Nth/2][k] << " " << Eth_tilde[NEW][i][Nth/2][k] << std::endl;
                //         }
                //         ofs_PML << std::endl;
                //     }
                // }

    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr - PML_L; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;

                // if(i == 1 && k == 50){
                //     ofs_Hph << n << " " << j << " " << Hphth[i][j][k] << " " << Hphr_tilde[NEW][i][j][k] << " " << Hphr[i][j][k] << " " << Hph[NEW][i][j][k] << " "
                //             << CHPHTH_00[j] << " " << CHPHTH_01[j] << " " << CHPHR_TILDE_00[i] << " " << CHPHR_TILDE_01[i] << " " << CHPHR_10[i] << " " << CHPHR_11[i] << " " << std::endl;
                // } 

            }
        }
    }

    ofs_Hph << std::endl;
    // #pragma omp parallel for collapse(3)
    for(int i = 1; i < Nr - PML_L; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Hth_tilde(double ***Hth_tilde, double ****Hth, double *CHTH_TILDE, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 0; i < PML_L + 1; i++){
    //     for(int j = 1; j <= Nth - 1; j++){
    //         for(int k = PML_L; k < Nph - PML_L; k++){
    //             Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L - 1; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L - 1; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    } 
}

void update_Hph_tilde(double ***Hph_tilde, double ****Hph, double *CHPH_TILDE, double **** check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 0; i < PML_L + 1; i++){
    //     for(int j = 0; j < Nth; j++){
    //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
    //             Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L - 1; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L - 1; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}
