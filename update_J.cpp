#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <string>

void update_Jr(double ****Jr, double ****Jth, double ****Jph, double ***Er, double ****Eth, double ****Eph, Eigen::Matrix3d ***S, Eigen::Matrix3d ***B, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    omp_set_num_threads(8);
    #pragma omp parallel for collapse(3)
    for(int i = Nr_iono_lower; i < Nr_iono_upper; i++){
        for(int j = Nth_iono_lower + 1; j <= Nth_iono_upper; j++){
            for(int k = Nph_iono_lower + 1; k <= Nph_iono_upper; k++){
                Jr[NEW][i][j][k] = S[i][j][k](0,0) * Jr[OLD][i][j][k]
                                + S[i][j][k](0,1) * (Jth[OLD][i][j][k] + Jth[OLD][i+1][j][k] + Jth[OLD][i][j-1][k] + Jth[OLD][i+1][j-1][k]) / 4.0
                                + S[i][j][k](0,2) * (Jph[OLD][i][j][k] + Jph[OLD][i+1][j][k] + Jph[OLD][i][j][k-1] + Jph[OLD][i][j+1][k-1]) / 4.0
                                + B[i][j][k](0,0) * Er[i][j][k]
                                + B[i][j][k](0,1) * (Eth[NEW][i][j][k] + Eth[NEW][i+1][j][k] + Eth[NEW][i][j-1][k] + Eth[NEW][i+1][j-1][k]) / 4.0
                                + B[i][j][k](0,2) * (Eph[NEW][i][j][k] + Eph[NEW][i+1][j][k] + Eph[NEW][i][j][k-1] + Eph[NEW][i][j+1][k-1]) / 4.0;
            }
        }
    }
}

void update_Jth(double ****Jr, double ****Jth, double ****Jph, double ***Er, double ****Eth, double ****Eph, Eigen::Matrix3d ***S, Eigen::Matrix3d ***B, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = Nr_iono_lower + 1; i <= Nr_iono_upper; i++){
        for(int j = Nth_iono_lower; j < Nth_iono_upper; j++){
            for(int k = Nph_iono_lower + 1; k <= Nph_iono_upper; k++){
                Jth[NEW][i][j][k] = S[i][j][k](1,0) * (Jr[OLD][i][j][k] + Jr[OLD][i][j+1][k] + Jr[OLD][i-1][j][k] + Jr[OLD][i-1][j+1][k]) / 4.0
                                + S[i][j][k](1,1) * Jth[OLD][i][j][k]
                                + S[i][j][k](1,2) * (Jph[OLD][i][j][k] + Jph[OLD][i][j+1][k] + Jph[OLD][i][j][k-1] + Jph[OLD][i][j+1][k-1]) / 4.0
                                + B[i][j][k](1,0) * (Er[i][j][k] + Er[i][j+1][k] + Er[i-1][j][k] + Er[i-1][j+1][k]) / 4.0
                                + B[i][j][k](1,1) * Eth[NEW][i][j][k]
                                + B[i][j][k](1,2) * (Eph[NEW][i][j][k] + Eph[NEW][i][j+1][k] + Eph[NEW][i][j][k-1] + Eph[NEW][i][j+1][k-1]) / 4.0;
            }
        }
    }
}

void update_Jph(double ****Jr, double ****Jth, double ****Jph, double ***Er, double ****Eth, double ****Eph, Eigen::Matrix3d ***S, Eigen::Matrix3d ***B, int n){
    int NEW = n % 2;
    int OLD = (n + 1) % 2;

    #pragma omp parallel for collapse(3)
    for(int i = Nr_iono_lower + 1; i <= Nr_iono_upper; i++){
        for(int j = Nth_iono_lower + 1; j <= Nth_iono_upper; j++){
            for(int k = Nph_iono_lower; k < Nph_iono_upper; k++){
                Jph[NEW][i][j][k] = S[i][j][k](2,0) * (Jr[OLD][i][j][k] + Jr[OLD][i][j][k+1] + Jr[OLD][i-1][j][k] + Jr[OLD][i-1][j][k+1]) / 4.0
                                + S[i][j][k](2,1) * (Jth[OLD][i][j][k] + Jth[OLD][i][j][k+1] + Jth[OLD][i][j-1][k] + Jth[OLD][i][j-1][k+1]) / 4.0
                                + S[i][j][k](2,2) * Jph[OLD][i][j][k]
                                + B[i][j][k](2,0) * (Er[i][j][k] + Er[i][j][k+1] + Er[i-1][j][k] + Er[i-1][j][k+1]) / 4.0
                                + B[i][j][k](2,1) * (Eth[NEW][i][j][k] + Eth[NEW][i][j][k+1] + Eth[NEW][i][j-1][k] + Eth[NEW][i][j-1][k+1]) / 4.0
                                + B[i][j][k](2,2) * Eph[NEW][i][j][k];
            }
        }
    }
}

