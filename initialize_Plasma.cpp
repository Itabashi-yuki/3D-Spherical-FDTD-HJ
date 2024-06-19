#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <string>

double cal_nu(double exp_nu){
    double nu = std::pow(10, exp_nu);
    return nu;
}

double cal_Ne(double exp_Ne){
    double Ne = std::pow(10, exp_Ne);
    return Ne;
}

double cal_omg_p(double Ne){
    return std::sqrt( Ne * CHARGE_e * CHARGE_e / ( MASS_e * EPS0 ) );
}

void initialize_Plasma(Eigen::Matrix3d ***S, Eigen::Matrix3d ***B){
    double exp_nu = 7.0;
    double exp_Ne = 8.0;
    double nu = cal_nu(exp_nu);
    double Ne = cal_Ne(exp_Ne);

    for(int i = Nr_iono_lower; i <= Nr_iono_upper; i++){
        for(int j = Nth_iono_lower; j <= Nth_iono_upper; j++){
            for(int k = Nph_iono_lower; k <= Nph_iono_upper; k++){
                double Omg_0 = 1.0 / dt + nu / 2.0;
                double Omg_0_prime = 1.0 / dt - nu / 2.0;
                double Omg_c = CHARGE_e * B0 / MASS_e;
                double Omg_p = cal_omg_p(Ne);

                Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                Eigen::Matrix3d A, R1, R2;
                Eigen::Vector3d b;

                b << std::sin(THETA) * std::cos(PHI), std::sin(THETA) * std::sin(PHI), std::cos(THETA);

                R2 << 0, b(2), -b(1),
                    -b(2), 0, b(0),
                    b(1), -b(0), 0;
                
                R1 = Omg_c / 2.0 * R2;

                A = Omg_0 * I + R1;

                S[i][j][k] = A.inverse() * ( Omg_0_prime * I - R1 );
                B[i][j][k] = EPS0 * Omg_p * Omg_p * A.inverse();
            
            }
        }
    }

}