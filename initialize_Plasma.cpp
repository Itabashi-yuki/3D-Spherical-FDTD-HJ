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

    std::string str_Hour = std::to_string(Hour);
    std::string str_Min = std::to_string(Min);
    str_Hour.insert(0,2 - str_Hour.length(), '0');
    str_Min.insert(0,2 - str_Min.length(), '0');

    std::ifstream ifs("./Ne/Ne_UT" + str_Hour + str_Min + ".dat");
    std::string line;
    std::getline(ifs, line); /* headerを読み飛ばす */
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
    // std::ofstream ofs_geomag( PATH + "data/" + global_dirName + "Geomag.dat",std::ios::app);
    for(int i = Nr_iono_lower; i <= Nr_iono_upper; i++){
        for(int j = Nth_iono_lower; j <= Nth_iono_upper; j++){
            for(int k = Nph_iono_lower; k <= Nph_iono_upper; k++){
                double alt = i * dr;
                double th_idx = j * dth;
                double ph_idx = k * dph;
                double F0, Dec, Inc;

                double nu = cal_nu(alt);
                // calc_geomagnetic_field(
                //     Year, Month, Day,
                //     z*1.0e-3, Tx_Latitude, Tx_Longitude,
                //     Dec, Inc, F0
                // ); 


                double Omg_0 = 1.0 / dt + nu / 2.0;
                double Omg_0_prime = 1.0 / dt - nu / 2.0;
                // double Omg_c = CHARGE_e * F0 / MASS_e;
                double Omg_c = CHARGE_e / MASS_e;
                double Omg_p = cal_omg_p(Ne);
                // double Omg_p = cal_omg_p(Ne_a[i]);

                Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                Eigen::Matrix3d A, R1, R2;
                Eigen::Vector3d b_car, b_sph;

                Eigen::Vector3d B0_components;
                
                // if(j==Nth_iono_lower+(Nr_iono/2) && k == Nph_iono_lower+(Nr_iono/2)){
                //     B0_components = transform_geomag(alt * 1.0e-3, th_idx, ph_idx, i, j, k);
                // }
                // B0_components = transform_geomag(alt * 1.0e-3, th_idx, ph_idx, i, j, k);
                
                // b << std::sin(THETA) * std::cos(PHI), std::sin(THETA) * std::sin(PHI), std::cos(THETA);

                // Eigen::Vector3d b;  // 3要素ベクトル
                b_car << std::sin(THETA) * std::cos(PHI), std::sin(THETA) * std::sin(PHI), std::cos(THETA);

                Eigen::Matrix3d Car2Sph;  // 3x3 行列
                Car2Sph << std::sin(th_idx)*std::cos(ph_idx), std::sin(th_idx)*std::sin(ph_idx), std::cos(th_idx),
                            std::cos(th_idx)*std::cos(ph_idx), std::cos(th_idx)*std::sin(ph_idx), -std::sin(th_idx),
                            -std::sin(ph_idx), std::cos(ph_idx), 0;

                b_sph = Car2Sph * b_car;  // 結果も3要素ベクトル

                R2 << 0, b_sph(2), -b_sph(1),
                        -b_sph(2), 0, b_sph(0),
                        b_sph(1), -b_sph(0), 0;

                // R2 << 0, B0_components(2), -B0_components(1),
                //     -B0_components(2), 0, B0_components(0),
                //     B0_components(1), -B0_components(0), 0;
                
                R1 = Omg_c / 2.0 * B0 * R2;

                A = Omg_0 * I + R1;

                S[i][j][k] = A.inverse() * ( Omg_0_prime * I - R1 );
                B[i][j][k] = EPS0 * Omg_p * Omg_p * A.inverse();
            
            }
        //    ofs_geomag << std::endl;
        }
    }

}