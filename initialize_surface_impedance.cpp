#include <iostream>
#include <cmath>

#include "fdtd3d.h"

double sig_s(double th, double ph){
    // return 1.0;
    // return 1.0e7;
    // return 1.0e-4;
    // if(ph <= 100){
    //     return 4.0;
    // } else if(ph <= 450){
    //     return 1.0e-4;
    // } else if(ph <= 600){
    //     return 4.0;
    // } else if(ph <= 700){
    //     return 1.0e-4;
    // } else {
    //     return 4.0;
    // }
    if(ph <= 100){
        return 1.0e-4;
    } else if(ph <= 450){
        return 4.0;
    } else if(ph <= 600){
        return 1.0e-4;
    } else if(ph <= 700){
        return 4.0;
    } else {
        return 1.0e-4;
    }
}

void intialize_surface_impedance(double **Rs, double **Ls){
    const double eps_r0 { Refractive_index(0.0) * Refractive_index(0.0) };

    for(int j = 0; j <= Nth; j++){
        for(int k = 0; k <= Nph; k++){
            double sig0 = sig_s(R0 * j *dth * 1e-3 - source_ph, R0 * k * dph * 1e-3);
            double r = std::sqrt(eps_r0 * eps_r0 + sig0 * sig0 / OMG / OMG / EPS0 / EPS0 );
            double th_arg = atan2(sig0, OMG * EPS0 * eps_r0);

            Rs[j][k] = Z0 / std::sqrt(r) * std::cos(th_arg / 2.0);
            Ls[j][k] = Z0 / OMG / std::sqrt(r) * std::sin(th_arg / 2.0);
        }
    }
}