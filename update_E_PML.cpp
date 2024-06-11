#include "fdtd3d.h"
#include <iostream>

void update_Er_PML(double ****Erth1, double ****Erth2, double ****Erph, double ***Er, double ***Hr, double ****Hth, double ****Hph,
                     double *CERTH1_00, double *CERTH1_01, double *CERPH_00, double *CERPH_01, double ****check, int n){

    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 1; k <= PML_L; k++){
                Erth1[NEW][i][j][k] = CERTH1_00[j] * Erth1[OLD][i][j][k] + CERTH1_01[j] / r(i + 0.5) / dth / EPS0 
                                    * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Erth2[NEW][i][j][k] = Erth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) / EPS0 
                                    * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Erph[NEW][i][j][k] = CERPH_00[k] * Erph[OLD][i][j][k] - CERPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph / EPS0 
                                    * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Er[i][j][k] = Erth1[NEW][i][j][k] + Erth2[NEW][i][j][k] + Erph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }


    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Erth1[NEW][i][j][k] = CERTH1_00[j] * Erth1[OLD][i][j][k] + CERTH1_01[j] / r(i + 0.5) / dth / EPS0 * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Erth2[NEW][i][j][k] = Erth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) / EPS0 * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Erph[NEW][i][j][k] = CERPH_00[k] * Erph[OLD][i][j][k] - CERPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph / EPS0 * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Er[i][j][k] = Erth1[NEW][i][j][k] + Erth2[NEW][i][j][k] + Erph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;        
            }
        }
    }

    for(int i = 0; i < PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Erth1[NEW][i][j][k] = CERTH1_00[j] * Erth1[OLD][i][j][k] + CERTH1_01[j] / r(i + 0.5) / dth / EPS0 * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Erth2[NEW][i][j][k] = Erth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) / EPS0 * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Erph[NEW][i][j][k] = CERPH_00[k] * Erph[OLD][i][j][k] - CERPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph / EPS0 * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Er[i][j][k] = Erth1[NEW][i][j][k] + Erth2[NEW][i][j][k] + Erph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;         
            }
        }
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Erth1[NEW][i][j][k] = CERTH1_00[j] * Erth1[OLD][i][j][k] + CERTH1_01[j] / r(i + 0.5) / dth / EPS0 * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Erth2[NEW][i][j][k] = Erth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) / EPS0 * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Erph[NEW][i][j][k] = CERPH_00[k] * Erph[OLD][i][j][k] - CERPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph / EPS0 * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Er[i][j][k] = Erth1[NEW][i][j][k] + Erth2[NEW][i][j][k] + Erph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Erth1[NEW][i][j][k] = CERTH1_00[j] * Erth1[OLD][i][j][k] + CERTH1_01[j] / r(i + 0.5) / dth / EPS0 * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Erth2[NEW][i][j][k] = Erth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) / EPS0 * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Erph[NEW][i][j][k] = CERPH_00[k] * Erph[OLD][i][j][k] - CERPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph / EPS0 * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Er[i][j][k] = Erth1[NEW][i][j][k] + Erth2[NEW][i][j][k] + Erph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L; i < Nr - PML_L; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Erth1[NEW][i][j][k] = CERTH1_00[j] * Erth1[OLD][i][j][k] + CERTH1_01[j] / r(i + 0.5) / dth / EPS0 * ( Hph[OLD][i][j][k] - Hph[OLD][i][j-1][k] );
                Erth2[NEW][i][j][k] = Erth2[OLD][i][j][k] + dt * cot(theta(j)) / 2.0 / r(i + 0.5) / EPS0 * ( Hph[OLD][i][j][k] + Hph[OLD][i][j-1][k] );
                Erph[NEW][i][j][k] = CERPH_00[k] * Erph[OLD][i][j][k] - CERPH_01[k] / r(i + 0.5) / std::sin(theta(j)) / dph / EPS0 * ( Hth[OLD][i][j][k] - Hth[OLD][i][j][k-1] );
                Er[i][j][k] = Erth1[NEW][i][j][k] + Erth2[NEW][i][j][k] + Erph[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }   
    // exit(0);
}

void update_Eth_PML(double ****Ethph, double ****Ethr, double ****Ethr_tilde, double ****Eth, double ***Hr, double ***Hph_tilde, double *CETHPH_00, double *CETHPH_01,
                     double *CETHR_10, double *CETHR_11, double *CETHR_TILDE_00, double *CETHR_TILDE_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n+1) % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Ethph[NEW][i][j][k] = CETHPH_00[k] * Ethph[OLD][i][j][k] + CETHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph / EPS0 * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Ethr_tilde[NEW][i][j][k] = CETHR_TILDE_00[i] * Ethr_tilde[OLD][i][j][k] - CETHR_TILDE_01[i] / dr / EPS0 * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Ethr[NEW][i][j][k] = CETHR_10[i] * Ethr[OLD][i][j][k] + CETHR_11[i] / dt * ( Ethr_tilde[NEW][i][j][k] - Ethr_tilde[OLD][i][j][k] );
                Eth[NEW][i][j][k] = Ethph[NEW][i][j][k] + Ethr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;

            }
        }
    }

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Ethph[NEW][i][j][k] = CETHPH_00[k] * Ethph[OLD][i][j][k] + CETHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph / EPS0 * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Ethr_tilde[NEW][i][j][k] = CETHR_TILDE_00[i] * Ethr_tilde[OLD][i][j][k] - CETHR_TILDE_01[i] / dr / EPS0 * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Ethr[NEW][i][j][k] = CETHR_10[i] * Ethr[OLD][i][j][k] + CETHR_11[i] / dt * ( Ethr_tilde[NEW][i][j][k] - Ethr_tilde[OLD][i][j][k] );
                Eth[NEW][i][j][k] = Ethph[NEW][i][j][k] + Ethr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= PML_L; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Ethph[NEW][i][j][k] = CETHPH_00[k] * Ethph[OLD][i][j][k] + CETHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph / EPS0 * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Ethr_tilde[NEW][i][j][k] = CETHR_TILDE_00[i] * Ethr_tilde[OLD][i][j][k] - CETHR_TILDE_01[i] / dr / EPS0 * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Ethr[NEW][i][j][k] = CETHR_10[i] * Ethr[OLD][i][j][k] + CETHR_11[i] / dt * ( Ethr_tilde[NEW][i][j][k] - Ethr_tilde[OLD][i][j][k] );
                Eth[NEW][i][j][k] = Ethph[NEW][i][j][k] + Ethr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Ethph[NEW][i][j][k] = CETHPH_00[k] * Ethph[OLD][i][j][k] + CETHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph / EPS0 * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Ethr_tilde[NEW][i][j][k] = CETHR_TILDE_00[i] * Ethr_tilde[OLD][i][j][k] - CETHR_TILDE_01[i] / dr / EPS0 * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Ethr[NEW][i][j][k] = CETHR_10[i] * Ethr[OLD][i][j][k] + CETHR_11[i] / dt * ( Ethr_tilde[NEW][i][j][k] - Ethr_tilde[OLD][i][j][k] );
                Eth[NEW][i][j][k] = Ethph[NEW][i][j][k] + Ethr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Ethph[NEW][i][j][k] = CETHPH_00[k] * Ethph[OLD][i][j][k] + CETHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph / EPS0 * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Ethr_tilde[NEW][i][j][k] = CETHR_TILDE_00[i] * Ethr_tilde[OLD][i][j][k] - CETHR_TILDE_01[i] / dr / EPS0 * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Ethr[NEW][i][j][k] = CETHR_10[i] * Ethr[OLD][i][j][k] + CETHR_11[i] / dt * ( Ethr_tilde[NEW][i][j][k] - Ethr_tilde[OLD][i][j][k] );
                Eth[NEW][i][j][k] = Ethph[NEW][i][j][k] + Ethr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Ethph[NEW][i][j][k] = CETHPH_00[k] * Ethph[OLD][i][j][k] + CETHPH_01[k] / r(i) / std::sin(theta(j+0.5)) / dph / EPS0 * ( Hr[i][j][k] - Hr[i][j][k-1] );
                Ethr_tilde[NEW][i][j][k] = CETHR_TILDE_00[i] * Ethr_tilde[OLD][i][j][k] - CETHR_TILDE_01[i] / dr / EPS0 * ( Hph_tilde[i][j][k] - Hph_tilde[i-1][j][k] );
                Ethr[NEW][i][j][k] = CETHR_10[i] * Ethr[OLD][i][j][k] + CETHR_11[i] / dt * ( Ethr_tilde[NEW][i][j][k] - Ethr_tilde[OLD][i][j][k] );
                Eth[NEW][i][j][k] = Ethph[NEW][i][j][k] + Ethr[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Eph_PML(double ****Ephr, double ****Ephr_tilde, double ****Ephth, double ****Eph, double ***Hr, double ***Hth_tilde,
                     double *CEPHR_10, double *CEPHR_11, double *CEPHR_TILDE_00, double *CEPHR_TILDE_01, double *CEPHTH_00, double *CEPHTH_01, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Ephth[NEW][i][j][k] = CEPHTH_00[j] * Ephth[OLD][i][j][k] - CEPHTH_01[j] / r(i) / dth / EPS0 * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Ephr_tilde[NEW][i][j][k] = CEPHR_TILDE_00[i] * Ephr_tilde[OLD][i][j][k] + CEPHR_TILDE_01[i] / dr / EPS0 * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Ephr[NEW][i][j][k] = CEPHR_10[i] * Ephr[OLD][i][j][k] + CEPHR_11[i] / dt * ( Ephr_tilde[NEW][i][j][k] - Ephr_tilde[OLD][i][j][k] );
                Eph[NEW][i][j][k] = Ephr[NEW][i][j][k] + Ephth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Ephth[NEW][i][j][k] = CEPHTH_00[j] * Ephth[OLD][i][j][k] - CEPHTH_01[j] / r(i) / dth / EPS0 * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Ephr_tilde[NEW][i][j][k] = CEPHR_TILDE_00[i] * Ephr_tilde[OLD][i][j][k] + CEPHR_TILDE_01[i] / dr / EPS0 * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Ephr[NEW][i][j][k] = CEPHR_10[i] * Ephr[OLD][i][j][k] + CEPHR_11[i] / dt * ( Ephr_tilde[NEW][i][j][k] - Ephr_tilde[OLD][i][j][k] );
                Eph[NEW][i][j][k] = Ephr[NEW][i][j][k] + Ephth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Ephth[NEW][i][j][k] = CEPHTH_00[j] * Ephth[OLD][i][j][k] - CEPHTH_01[j] / r(i) / dth / EPS0 * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Ephr_tilde[NEW][i][j][k] = CEPHR_TILDE_00[i] * Ephr_tilde[OLD][i][j][k] + CEPHR_TILDE_01[i] / dr / EPS0 * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Ephr[NEW][i][j][k] = CEPHR_10[i] * Ephr[OLD][i][j][k] + CEPHR_11[i] / dt * ( Ephr_tilde[NEW][i][j][k] - Ephr_tilde[OLD][i][j][k] );
                Eph[NEW][i][j][k] = Ephr[NEW][i][j][k] + Ephth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Ephth[NEW][i][j][k] = CEPHTH_00[j] * Ephth[OLD][i][j][k] - CEPHTH_01[j] / r(i) / dth / EPS0 * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Ephr_tilde[NEW][i][j][k] = CEPHR_TILDE_00[i] * Ephr_tilde[OLD][i][j][k] + CEPHR_TILDE_01[i] / dr / EPS0 * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Ephr[NEW][i][j][k] = CEPHR_10[i] * Ephr[OLD][i][j][k] + CEPHR_11[i] / dt * ( Ephr_tilde[NEW][i][j][k] - Ephr_tilde[OLD][i][j][k] );
                Eph[NEW][i][j][k] = Ephr[NEW][i][j][k] + Ephth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Ephth[NEW][i][j][k] = CEPHTH_00[j] * Ephth[OLD][i][j][k] - CEPHTH_01[j] / r(i) / dth / EPS0 * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Ephr_tilde[NEW][i][j][k] = CEPHR_TILDE_00[i] * Ephr_tilde[OLD][i][j][k] + CEPHR_TILDE_01[i] / dr / EPS0 * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Ephr[NEW][i][j][k] = CEPHR_10[i] * Ephr[OLD][i][j][k] + CEPHR_11[i] / dt * ( Ephr_tilde[NEW][i][j][k] - Ephr_tilde[OLD][i][j][k] );
                Eph[NEW][i][j][k] = Ephr[NEW][i][j][k] + Ephth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Ephth[NEW][i][j][k] = CEPHTH_00[j] * Ephth[OLD][i][j][k] - CEPHTH_01[j] / r(i) / dth / EPS0 * ( Hr[i][j][k] - Hr[i][j-1][k] );
                Ephr_tilde[NEW][i][j][k] = CEPHR_TILDE_00[i] * Ephr_tilde[OLD][i][j][k] + CEPHR_TILDE_01[i] / dr / EPS0 * ( Hth_tilde[i][j][k] - Hth_tilde[i-1][j][k] );
                Ephr[NEW][i][j][k] = CEPHR_10[i] * Ephr[OLD][i][j][k] + CEPHR_11[i] / dt * ( Ephr_tilde[NEW][i][j][k] - Ephr_tilde[OLD][i][j][k] );
                Eph[NEW][i][j][k] = Ephr[NEW][i][j][k] + Ephth[NEW][i][j][k];
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }    
}

void update_Eth_tilde(double ****Eth_tilde, double ****Eth, double *CETH_TILDE_00, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= PML_L; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Eph_tilde(double ****Eph_tilde, double ****Eph, double *CEPH_TILDE_00, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = 1; i <= PML_L; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    for(int i = PML_L + 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L + 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

}