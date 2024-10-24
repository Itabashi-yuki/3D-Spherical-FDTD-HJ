#include "fdtd3d.h"

double source_J(double t){
    return  ( t - t0 ) / sigma_J * std::exp( -(t - t0) * (t - t0) / 2.0 / sigma_J / sigma_J );

    // double Jr = 0.0;
    // if(0 <= t && t < t0){
    //     Jr = sin( 2 * M_PI * f0 * t) * std::exp( -( t - t0) * ( t - t0) / sigma_J / sigma_J );
    // } else if (t0 < t){
    //     Jr = sin(2 * M_PI * f0 * t);
    // }
    // return Jr;
    
    // double Jr = std::exp(-((t - t0) * (t - t0)) / (2.0 * sigma_J * sigma_J));
    // return Jr * std::sin(OMG * t);

    // return 0;
}