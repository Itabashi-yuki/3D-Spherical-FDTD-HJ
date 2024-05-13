#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <string>


void output_E(double ***Er, double ****Eth, double ****Eph, double ***Hr, double ****Hth, double ****Hph, int n, int n0){
    int NEW = n % 2;
    if(n == 1 || n % n0 == 0 || n == Nt - 1){
        std::ofstream ofs_E("./data/" + global_dirName + "/E_" + std::to_string(n) + ".dat");

        ofs_E << "#" << n0 * dt << "[s]ごと出力, 以下は" << n*dt << "[s]における結果" <<  std::endl;

        // for(int i = 0; i < Nr; i+=5){
        //     for(int j = 0; j < Nth; j+=5){
        //         for(int k = 0; k < Nph; k+=5){
        //             // ofs_Ez << i * dx * 1e-3 << " " << j * dy * 1e-3 << " " << Ez[i][j] << std::endl;
        //             // ofs_E << (R0 + i * dr) * 1e-3 << " " << j * R0 * dth * 1e-3 << " " << k * R0 * dph * 1e-3 << " " 
        //             // << Er[i][j][k] << " " << Eth[i][j][k] << " " << Eph[i][j][k] << " " 
        //             // << std::sqrt(Er[i][j][k]*Er[i][j][k] + Eth[i][j][k]*Eth[i][j][k] + Eph[i][j][k]*Eph[i][j][k]) << std::endl;
        //             ofs_E << j * R0 * dth * 1e-3 << " " << Er[i][j][k] << " " << Eth[i][j][k] << " " << Eph[i][j][k] << " " 
        //             << Hr[i][j][k] << " " << Hth[i][j][k] << " " << Hph[i][j][k] << std::endl;
        //         }
        //     ofs_E << std::endl;
        //     }
        // }

            for(int k = 0; k < Nph; k+=2){
                // ofs_E << j * R0 * dth * 1e-3 << " " << Er[Nr / 2][j][Nph / 2] << std::endl;
                ofs_E << k * R0 * dth * 1e-3 << " " << Er[Nr / 2][Nth / 2][k] << " " << Eth[NEW][Nr / 2][Nth / 2][k] << " " << Eph[NEW][Nr / 2][Nth / 2][k] << " " << Hr[Nr / 2][Nth / 2][k] 
                << " " << Hth[NEW][Nr / 2][Nth / 2][k] << " " << Hph[NEW][Nr / 2][Nth / 2][k] << std::endl;
                // ofs_E << k * R0 * dph * 1e-3 << " " << Er[Nr / 2][Nth / 2][k] << std::endl;
            }

    ofs_E.close();

    }
}