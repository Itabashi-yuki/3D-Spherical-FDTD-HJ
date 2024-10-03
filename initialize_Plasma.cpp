#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <string>
#include <geomag.h>
#include <Vector3d.h>
#include <cmath>

double cal_nu(double z){
    return 4.303e11 * exp( -0.1622 * z * 1.0e-3 );
}

double cal_Ne(double exp_Ne){
    double Ne = std::pow(10, exp_Ne);
    return 0.0;
    // return Ne;
}

double cal_omg_p(double Ne){
    return std::sqrt( Ne * CHARGE_e * CHARGE_e / ( MASS_e * EPS0 ) );
}

void initialize_Plasma(Eigen::Matrix3d ***S, Eigen::Matrix3d ***B){
    // double exp_nu = 7.0;
    // double exp_Ne = 9.0;
    // double nu = cal_nu(exp_nu);
    double Ne = cal_Ne(exp_Ne);
    double *Ne_a = allocate_1d(Nr, 0.0);
    double *Alt = allocate_1d(Nr, 0.0);
    // std::ifstream ifs("interpolate_Ne_UT1700.dat");
    // std::ifstream ifs("interpolated_Ne_UT1700_firi_05km.dat");
    std::ifstream ifs("interpolate_Ne_UT1700.dat");
    for(int i = Nr_iono_lower; i <= Nr_iono_upper; i++){
        ifs >> Alt[i] >> Ne_a[i];
    }
    std::ofstream ofs( PATH + "data/" + global_dirName + "Ne.dat");
    for(int i = Nr_iono_lower; i <= Nr_iono_upper; i++){
        ofs << i * dr * 1e-3 << " " << Ne_a[i] << std::endl;
    }
    // exit(0);
    // double Ne_max_exp = 12.0;
    // double Ne_min_exp = 8.0;

    // for(int i = Nr_iono_lower; i <= Nr_iono_upper; i++){
    //     int r = i - Nr_iono_lower;
    //     Ne_a[i] = std::pow(10, (Ne_max_exp * r + Ne_min_exp * (Nr_iono - r)) / Nr_iono);
    //     ofs << i * dr * 1e-3 << " " << Ne_a[i] << std::endl;
    // }

    // exit(0);

    for(int i = Nr_iono_lower; i <= Nr_iono_upper; i++){
        for(int j = Nth_iono_lower; j <= Nth_iono_upper; j++){
            for(int k = Nph_iono_lower; k <= Nph_iono_upper; k++){
                double z = i * dr;
                double th = j * dth;
                double ph = k * dph;
                double F0, Dec, Inc;

                double nu = cal_nu(z);
                calc_geomagnetic_field(
                    Year, Month, Day,
                    z*1.0e-3, Tx_Latitude, Tx_Longitude,
                    Dec, Inc, F0
                ); 

                double Omg_0 = 1.0 / dt + nu / 2.0;
                double Omg_0_prime = 1.0 / dt - nu / 2.0;
                double Omg_c = CHARGE_e * F0 / MASS_e;
                // double Omg_p = cal_omg_p(Ne);
                double Omg_p = cal_omg_p(Ne_a[i]);

                Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                Eigen::Matrix3d A, R1, R2;
                Eigen::Vector3d b;

                Eigen::Vector3d B0_components;
                B0_components = transform_geomag(Inc, Dec, th, ph);
                // b << std::sin(THETA) * std::cos(PHI), std::sin(THETA) * std::sin(PHI), std::cos(THETA);

                R2 << 0, B0_components(2), -B0_components(1),
                    -B0_components(2), 0, B0_components(0),
                    B0_components(1), -B0_components(0), 0;
                
                R1 = Omg_c / 2.0 * R2;

                A = Omg_0 * I + R1;

                S[i][j][k] = A.inverse() * ( Omg_0_prime * I - R1 );
                B[i][j][k] = EPS0 * Omg_p * Omg_p * A.inverse();
            
            }
        }
    }

}