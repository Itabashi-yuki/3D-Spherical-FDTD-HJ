#include "fdtd3d.h"
#include <iostream>
#include <cmath>

void update_Hr_PML(double ***Hr, double ***Hrth1, double ***Hrth2, double ***Hrph,
                     double ****Eth, double ****Eph, double *CHRTH1_00, double *CHRTH1_01,
                      double *CHRPH_00, double *CHRPH_01, double ****check, int n){
    int NEW = n % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 0; k < PML_L; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Hrth1[i][j][k] = CHRTH1_00[j] * Hrth1[i][j][k] - CHRTH1_01[j] / MU0 / r(i) / dth * (Eph[NEW][i][j+1][k] - Eph[NEW][i][j][k]);
                Hrth2[i][j][k] = Hrth2[i][j][k] - dt * cot(theta(j+0.5)) / 2.0 / MU0 / r(i) * (Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k]);
                Hrph[i][j][k] = CHRPH_00[k] * Hrph[i][j][k] + CHRPH_01[k] / MU0 / r(i) / std::sin(theta(j + 0.5)) / dph * (Eth[NEW][i][j][k+1] - Eth[NEW][i][j][k]);
                Hr[i][j][k] = Hrth1[i][j][k] + Hrth2[i][j][k] + Hrph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= PML_L; i++){
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

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
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

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
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

void update_Hth_PML(double ****Hth, double ***Hthr, double ***Hthph, double ****Hthr_tilde, double ***Er, double ****Eph_tilde, double *CHTHPH_00,
                 double *CHTHPH_01, double *CHTHR_TILDE_00, double *CHTHR_TILDE_01, double *CHTHR_10, double *CHTHR_11, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hthph[i][j][k] = CHTHPH_00[k] * Hthph[i][j][k] - CHTHPH_01[k] / MU0 / r(i+0.5) / std::sin(theta(j)) / dph * (Er[i][j][k+1] - Er[i][j][k]);
                Hthr_tilde[NEW][i][j][k] = CHTHR_TILDE_00[i] * Hthr_tilde[OLD][i][j][k] + CHTHR_TILDE_01[i] / MU0 / dr * (Eph_tilde[NEW][i+1][j][k] - Eph_tilde[NEW][i][j][k]);
                Hthr[i][j][k] = CHTHR_10[i] * Hthr[i][j][k] + CHTHR_11[i] / dt * (Hthr_tilde[NEW][i][j][k] - Hthr_tilde[OLD][i][j][k]);
                Hth[NEW][i][j][k] = Hthr[i][j][k] + Hthph[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
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

void update_Hph_PML(double ****Hph, double ***Hphr, double ***Hphth, double ****Hphr_tilde, double ***Er, double ****Eth_tilde,
                     double *CHPHTH_00, double *CHPHTH_01, double *CHPHR_TILDE_00, double *CHPHR_TILDE_01, double *CHPHR_10, double *CHPHR_11, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 0; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < PML_L; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hphth[i][j][k] = CHPHTH_00[j] * Hphth[i][j][k] + CHPHTH_01[j] / MU0 / r(i + 0.5) / dth * (Er[i][j+1][k] - Er[i][j][k]);
                Hphr_tilde[NEW][i][j][k] = CHPHR_TILDE_00[i] * Hphr_tilde[OLD][i][j][k] - CHPHR_TILDE_01[i] / MU0 / dr * (Eth_tilde[NEW][i+1][j][k] - Eth_tilde[NEW][i][j][k]);
                Hphr[i][j][k] = CHPHR_10[i] * Hphr[i][j][k] + CHPHR_11[i] / dt * (Hphr_tilde[NEW][i][j][k] - Hphr_tilde[OLD][i][j][k]);
                Hph[NEW][i][j][k] = Hphth[i][j][k] + Hphr[i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
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

    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Hth_tilde[i][j][k] = Hth_tilde[i][j][k] + r(i + 0.5) * ( Hth[NEW][i][j][k] - Hth[OLD][i][j][k] ) + CHTH_TILDE[i] * ( Hth[NEW][i][j][k] + Hth[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
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

    for(int i = 0; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 0; i < PML_L; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Hph_tilde[i][j][k] = Hph_tilde[i][j][k] + r(i+0.5) * ( Hph[NEW][i][j][k] - Hph[OLD][i][j][k] ) + CHPH_TILDE[i] * ( Hph[NEW][i][j][k] + Hph[OLD][i][j][k] );
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}