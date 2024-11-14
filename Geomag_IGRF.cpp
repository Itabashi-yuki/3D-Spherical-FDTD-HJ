#include <iostream>
#include <fstream>
#include <cmath>

#include <Vector3d.h>
#include <geomag.h>
#include <eigen3/Eigen/Dense>
#include "fdtd3d.h"

Eigen::Vector3d transform_geomag(double alt, double th_idx, double ph_idx, int i, int j, int k){

    /* 送信点ベクトルnR_T, 受信点ベクトルnR_R */
    AndoLab::Vector3d <double> nR_T(1.0, (90.0 - Tx_Latitude)*AndoLab::DEG2RAD, Tx_Longitude*AndoLab::DEG2RAD, AndoLab::coordinate::Spherical);
    AndoLab::Vector3d <double> nR_R(1.0, (90.0 - Rx_Latitude)*AndoLab::DEG2RAD, Rx_Longitude*AndoLab::DEG2RAD, AndoLab::coordinate::Spherical);

    /* シミュレーション座標の基本単位ベクトル x, y, z */
    AndoLab::Vector3d <double> x = nR_T;
    AndoLab::Vector3d <double> z = (nR_T * nR_R).n();
    AndoLab::Vector3d <double> y = z * x;

    /* 送信点方向と受信点方向のなす角 */
    double ph_r = angle_between(nR_T, nR_R);
    // std::cout << ph_r * AndoLab::RAD2DEG << std::endl;
    // exit(0);
    /* 解析領域中央の点 */
    // double th { M_PI / 2.0 };
    // double ph { ph_r / 2.0 };
    
    double th = M_PI / 2.0 - thR + th_idx;
    double ph = ph_idx;
    // std::cout << th << " " << ph << std::endl;

    AndoLab::Vector3d <double> r_v = 
    std::sin(th) * std::cos(ph) * x + std::sin(th) * std::sin(ph) * y + std::cos(th) * z;
    
    AndoLab::Vector3d <double> th_v = 
    std::cos(th) * std::cos(ph) * x + std::cos(th) * std::sin(ph) * y - std::sin(th) * z;
    
    AndoLab::Vector3d <double> ph_v = 
    -std::sin(ph) * x + std::cos(ph) * y;

    // std::cout << r_v.x() << " " << r_v.y() << " " << r_v.z() << std::endl;
    double r_X = r_v.x();
    double r_Y = r_v.y();
    double r_Z = r_v.z();

    /* 任意の点における緯度経度 */
    double Th = std::atan2(std::sqrt(r_X * r_X + r_Y * r_Y), r_Z);
    double Ph = std::atan2(r_Y, r_X);
    // std::cout << (M_PI / 2.0 - Th) * AndoLab::RAD2DEG << " " << Ph * AndoLab::RAD2DEG << std::endl;
    // exit(0);

    double Inc, Dec, F0;
    calc_geomagnetic_field(
        Year, Month, Day,
        alt, (M_PI / 2.0 - Th)*AndoLab::RAD2DEG, Ph*AndoLab::RAD2DEG,
        Dec, Inc, F0
    ); 

    std::ofstream ofs_geomag( PATH + "data/" + global_dirName + "Geomag.dat",std::ios::app);

    // if((i-12) % 10 == 0){
    //     ofs_geomag << alt-12 << " " << j << " " << k << " " << (M_PI/2.0 - Th)*AndoLab::RAD2DEG << " " << Ph*AndoLab::RAD2DEG 
    //                 << " " << Dec << " " << Inc << " "<< F0 << " " << std::sin((90.0 - Dec)*AndoLab::DEG2RAD) << " " 
    //                 << std::cos((90.0 - Dec)*AndoLab::DEG2RAD) << " " << std::sin(-Inc*AndoLab::DEG2RAD) <<std::endl; 
    // }
    // ofs_geomag << (M_PI/2.0 -Th)*AndoLab::RAD2DEG << " " << Ph*AndoLab::RAD2DEG << " " << Dec << " " << Inc << " " 
    //             << std::sin(M_PI/2.0 + Ph) << " " << std::cos(M_PI/2.0 + Ph) << std::endl;
        // exit(0);

    /* 地理座標における基本単位ベクトル R_v, Th_v, Ph_v
        いずれも r_v, th_v, ph_vを用いて表わせる */
    AndoLab::Vector3d <double> R_v(1.0, Th, Ph, AndoLab::coordinate::Spherical);
    AndoLab::Vector3d <double> Th_v = R_v.theta_vector();
    AndoLab::Vector3d <double> Ph_v = R_v.phi_vector();

    /* 地磁気を求める */
    AndoLab::Vector3d <double> B0
    = F0 * (-std::sin(Inc) * R_v - std::cos(Inc) * std::cos(Dec) * Th_v + std::cos(Inc) * std::sin(Dec) * Ph_v);
    // std::cout << B0.abs() << std::endl;
    
    // r, th, ph成分を取り出して行ベクトルに保存
    Eigen::Vector3d B_components;
    B_components << B0 % r_v, B0 % th_v, B0 % ph_v;  // それぞれの成分を行ベクトルに格納

    // B0.convert2cartesian();
    // Eigen::Vector3d B_components;
    // B_components << B0 % x, B0 % y, B0 % z;  // それぞれの成分を行ベクトルに格納

    return B_components;  // 地磁気方向ベクトルを返す
}