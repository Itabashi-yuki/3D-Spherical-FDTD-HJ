#include "fdtd3d.h"

double Jr(double t){
    return  - (t - t0) / sigma_J * std::exp( -(t - t0) * (t - t0) / 2.0 / sigma_J / sigma_J );

}