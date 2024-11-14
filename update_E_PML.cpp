#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <string>

void update_Er_PML(double ***Er, double ****Dr, double ****Hth, double ****Hph, double ****Jr, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;  

    // std::ofstream ofs(PATH + "data/" + global_dirName + "Er_check.dat",std::ios::app);
    omp_set_num_threads(10);
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 1; k <= PML_L; k++){
                Er[i][j][k] = Er[i][j][k] + 1.0 / EPS0 * ( Dr[NEW][i][j][k] - Dr[OLD][i][j][k] ) - dt / EPS0 * Jr[OLD][i][j][k]; 
                // if(i == 0){
			// ofs << n << " " << j << " " << k << " "  << Er[i][j][k] << " " <<  Dr[NEW][i][j][k] << " " << Dr[OLD][i][j][k] << " " << Dr[NEW][i][j][k] - Dr[OLD][i][j][k] << std::endl;
		// }]
        // check[NEW][i][j][k] += 1.0;
                // Er[PML_L][j][k] = 0.0;
            }
        }
    }
	// ofs << std::endl;

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Er[i][j][k] = Er[i][j][k] + 1.0 / EPS0 * ( Dr[NEW][i][j][k] - Dr[OLD][i][j][k] ) - dt / EPS0 * Jr[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
                // Er[PML_L][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 0; i < PML_L; i++){
    //     for(int j = 1; j <= Nth - 1; j++){
    //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
    //             Er[i][j][k] = Er[i][j][k] + 1.0 / EPS0 * ( Dr[NEW][i][j][k] - Dr[OLD][i][j][k] ) - dt / EPS0 * Jr[OLD][i][j][k]; 
    //             // check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i < Nr; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Er[i][j][k] = Er[i][j][k] + 1.0 / EPS0 * ( Dr[NEW][i][j][k] - Dr[OLD][i][j][k] ) - dt / EPS0 * Jr[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }


    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Er[i][j][k] = Er[i][j][k] + 1.0 / EPS0 * ( Dr[NEW][i][j][k] - Dr[OLD][i][j][k] ) - dt / EPS0 * Jr[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nr - PML_L; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Er[i][j][k] = Er[i][j][k] + 1.0 / EPS0 * ( Dr[NEW][i][j][k] - Dr[OLD][i][j][k] ) - dt / EPS0 * Jr[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Eth_PML(double ****Eth, double ****Dth, double ***Hr, double ****Hph, double ****Jth, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;
    
    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= PML_L; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + 1.0 / EPS0 * ( Dth[NEW][i][j][k] - Dth[OLD][i][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
                // Eth[NEW][PML_L + 1][j][k] = 0.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = Nph - PML_L; k <= Nph - 1; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + 1.0 / EPS0 * ( Dth[NEW][i][j][k] - Dth[OLD][i][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
                // Eth[NEW][PML_L + 1][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i <= PML_L; i++){
    //     for(int j = 0; j < Nth; j++){
    //         for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
    //             Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + 1.0 / EPS0 * ( Dth[NEW][i][j][k] - Dth[OLD][i][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k]; 
    //             check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + 1.0 / EPS0 * ( Dth[NEW][i][j][k] - Dth[OLD][i][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }


    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = 0; j < PML_L; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + 1.0 / EPS0 * ( Dth[NEW][i][j][k] - Dth[OLD][i][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j < Nth; j++){
            for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
                Eth[NEW][i][j][k] = Eth[OLD][i][j][k] + 1.0 / EPS0 * ( Dth[NEW][i][j][k] - Dth[OLD][i][j][k] ) - dt / EPS0 * Jth[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Eph_PML(double ****Eph, double ****Dph, double ***Hr, double ****Hth, double ****Jph, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < PML_L; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + 1.0 / EPS0 * ( Dph[NEW][i][j][k] - Dph[OLD][i][j][k] ) - dt / EPS0 * Jph[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
                // Eph[NEW][PML_L + 1][j][k] = 0.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = Nph - PML_L; k < Nph; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + 1.0 / EPS0 * ( Dph[NEW][i][j][k] - Dph[OLD][i][j][k] ) - dt / EPS0 * Jph[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
                // Eph[NEW][PML_L + 1][j][k] = 0.0;
            }
        }
    }

    // #pragma omp parallel for collapse(3)
    // for(int i = 1; i <= PML_L; i++){
    //     for(int j = 1; j <= Nth - 1; j++){
    //         for(int k = PML_L; k < Nph - PML_L; k++){
    //             Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + 1.0 / EPS0 * ( Dph[NEW][i][j][k] - Dph[OLD][i][j][k] ) - dt / EPS0 * Jph[OLD][i][j][k]; 
    //             check[NEW][i][j][k] += 1.0;
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(3)
    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + 1.0 / EPS0 * ( Dph[NEW][i][j][k] - Dph[OLD][i][j][k] ) - dt / EPS0 * Jph[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }


    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = 1; j <= PML_L; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + 1.0 / EPS0 * ( Dph[NEW][i][j][k] - Dph[OLD][i][j][k] ) - dt / EPS0 * Jph[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - PML_L - 1; i++){
        for(int j = Nth - PML_L; j <= Nth - 1; j++){
            for(int k = PML_L; k < Nph - PML_L; k++){
                Eph[NEW][i][j][k] = Eph[OLD][i][j][k] + 1.0 / EPS0 * ( Dph[NEW][i][j][k] - Dph[OLD][i][j][k] ) - dt / EPS0 * Jph[OLD][i][j][k]; 
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

void update_Eth_tilde(double ****Eth_tilde, double ****Eth, double *CETH_TILDE_00, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 0; j < Nth; j++){
            for(int k = 1; k <= Nph - 1; k++){
                Eth_tilde[NEW][i][j][k] = Eth_tilde[OLD][i][j][k] + r(i) * (Eth[NEW][i][j][k] - Eth[OLD][i][j][k]) 
                                        + CETH_TILDE_00[i] * (Eth[NEW][i][j][k] + Eth[OLD][i][j][k]);
                // check[NEW][i][js][k] += 1.0;
            }
        }
    }
}

void update_Eph_tilde(double ****Eph_tilde, double ****Eph, double *CEPH_TILDE_00, double ****check, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = 1; i <= Nr - 1; i++){
        for(int j = 1; j <= Nth - 1; j++){
            for(int k = 0; k < Nph; k++){
                Eph_tilde[NEW][i][j][k] = Eph_tilde[OLD][i][j][k] + r(i) * (Eph[NEW][i][j][k] - Eph[OLD][i][j][k]) 
                                        + CEPH_TILDE_00[i] * (Eph[NEW][i][j][k] + Eph[OLD][i][j][k]);
                // check[NEW][i][j][k] += 1.0;
            }
        }
    }
}

