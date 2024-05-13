#include "fdtd3d.h"

int cal_obs_n0(){
    double n0 = 0.0;
    for(int n = 1; n < Nt; n++){
        if(int(n*dt / obs_t_step) == 1 ){
            n0 = n;
            break;
        }
    }
    return n0;
}