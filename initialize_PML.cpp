#include "fdtd3d.h"
#include <iostream>
#include <fstream>

double PML_sigma_r(double i){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / dr * std::log(PML_R);
    return sigma_max * std::pow((PML_L - i) * dr / PML_L / dr, PML_M);
    // return sigma_max * std::pow((PML_L * dr - r(i)) / PML_L / dr, PML_M);
}

double PML_sigma_th(double j){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / dth * std::log(PML_R);
    return sigma_max * std::pow((PML_L - j) * dth / PML_L / dth, PML_M);
}

double PML_sigma_ph(double k){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / dph * std::log(PML_R);
    return sigma_max * std::pow((PML_L - k) * dph / PML_L / dph, PML_M);
}

double PML_sigma_r_tilde(double i){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / dr * std::log(PML_R);
    // return sigma_max * ( PML_L * dr / (PML_M + 1.0)) * std::pow(( (PML_L - i) / PML_L ), (PML_M + 1));
    return sigma_max * (  PML_L  / (PML_M + 1.0)) * std::pow(( (PML_L - i) * dr / PML_L / dr ), (PML_M + 1.0));
    // return sigma_max * ( - PML_L / ( PML_M + 1.0) ) * std::pow(  ( PML_L * dr - r(i) ) / ( PML_L * dr ), (PML_M + 1.0));
}

double PML_C00(double sigma){
    return ( 1.0 / dt - sigma / 2.0 ) / ( 1.0 / dt + sigma / 2.0); 
}

double PML_C01(double sigma){
    return 1.0 / ( 1.0 / dt + sigma / 2.0); 
}

double PML_C10(double r, double sigma){
    return ( r / dt - sigma / 2.0) / ( r / dt + sigma / 2.0 );
}

double PML_C11(double r, double sigma){
    return 1.0 / ( r / dt + sigma / 2.0 );
}

void initialize_PML(double *CERTH1_00, double *CERTH1_01, double *CERPH_00, double *CERPH_01, double *CETHPH_00, double *CETHPH_01, double *CETHR_10, double *CETHR_11,
                     double *CETHR_TILDE_00, double *CETHR_TILDE_01, double *CEPHR_10, double *CEPHR_11, double *CEPHR_TILDE_00, double *CEPHR_TILDE_01, double *CEPHTH_00,
                     double *CEPHTH_01, double *CETH_TILDE_00, double *CETH_TILDE_01,
                     double *CEPH_TILDE_00, double *CEPH_TILDE_01, double *CHRTH1_00, double *CHRTH1_01, double *CHRPH_00, double *CHRPH_01, double *CHTHPH_00,
                     double *CHTHPH_01, double *CHTHR_TILDE_00, double *CHTHR_TILDE_01, double *CHTHR_10, double *CHTHR_11,
                     double *CHPHTH_00, double *CHPHTH_01, double *CHPHR_TILDE_00, double *CHPHR_TILDE_01, double *CHPHR_10, double *CHPHR_11, double *CHTH_TILDE, double *CHPH_TILDE){
    double *sigma_r_tilde_E = allocate_1d(Nr, 0.0);        
    double *sigma_th_E = allocate_1d(Nth, 0.0);
    double *sigma_ph_E = allocate_1d(Nph, 0.0);
    double *sigma_r_E = allocate_1d(Nr, 0.0);
    double *sigma_r_H = allocate_1d(Nr, 0.0);
    double *sigma_th_H = allocate_1d(Nth, 0.0);
    double *sigma_ph_H = allocate_1d(Nph, 0.0);
    double *sigma_r_tilde_H = allocate_1d(Nr, 0.0);

    for(int i = 1; i <= PML_L; i++){
        sigma_r_tilde_E[i] = PML_sigma_r_tilde(i);
    }

    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        sigma_r_tilde_E[i] = PML_sigma_r_tilde(Nr - i);
    }

    for(int i = 1; i <= PML_L; i++){
        sigma_r_E[i] = PML_sigma_r(i);
    }

    for(int i = Nr - PML_L; i <= Nr - 1; i++){
        sigma_r_E[i] = PML_sigma_r(Nr - i);
    }

    for(int j = 1; j <= PML_L; j++){
        sigma_th_E[j] = PML_sigma_th(j);
    }

    for(int j = Nth - PML_L; j <= Nth - 1; j++){
        sigma_th_E[j] = PML_sigma_th(Nth - j);
    }

    for(int k = 1; k <= PML_L; k++){
        sigma_ph_E[k] = PML_sigma_ph(k);
    }

    for(int k = Nph - PML_L; k <= Nph - 1; k++){
        sigma_ph_E[k] = PML_sigma_ph(Nph - k);
    }

    for(int i = 0; i < PML_L; i++){
        sigma_r_H[i] = PML_sigma_r(i + 0.5);
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        sigma_r_H[i] = PML_sigma_r(Nr - (i + 0.5));
    }

    for(int j = 0; j < PML_L; j++){
        sigma_th_H[j] = PML_sigma_th(j + 0.5);
    }

    for(int j = Nth - PML_L; j < Nth; j++){
        sigma_th_H[j] = PML_sigma_th(Nth - (j + 0.5));
    }
    
    for(int k = 0; k < PML_L; k++){
        sigma_ph_H[k] = PML_sigma_ph(k + 0.5);
    }

    for(int k = Nph - PML_L; k < Nph; k++){
        sigma_ph_H[k] = PML_sigma_ph(Nph - (k + 0.5));
    }

    for(int i = 0; i < PML_L; i++){
        sigma_r_tilde_H[i] = PML_sigma_r_tilde(i + 0.5);
    }

    for(int i = Nr - PML_L; i < Nr; i++){
        sigma_r_tilde_H[i] = PML_sigma_r_tilde(Nr - (i + 0.5));
    }

    std::ofstream ofs_sigma("./data/" + global_dirName + "/sigma.dat");
    for(int j = 0; j < Nth; j++){
        ofs_sigma << j << " " << sigma_r_E[j] << " " << sigma_th_E[j] << " " << sigma_ph_E[j] << " " << sigma_r_H[j] << " " << sigma_th_H[j] << " " 
                    << sigma_ph_H[j] << " " << sigma_r_tilde_E[j] << " " << sigma_r_tilde_H[j] << std::endl;
    }

    ofs_sigma.close();
    // exit(0);

    for(int j = 1; j <= Nth - 1; j++){
        CERTH1_00[j] = PML_C00(sigma_th_E[j]);
        CERTH1_01[j] = PML_C01(sigma_th_E[j]);
    }

    for(int k = 1; k <= Nph - 1; k++){
        CERPH_00[k] = PML_C00(sigma_ph_E[k]);
        CERPH_01[k] = PML_C01(sigma_ph_E[k]);
    }

    for(int k = 1; k <= Nph - 1; k++){
        CETHPH_00[k] = PML_C00(sigma_ph_E[k]);
        CETHPH_01[k] = PML_C01(sigma_ph_E[k]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CETHR_10[i] = PML_C10(r(i), sigma_r_tilde_E[i]);
        CETHR_11[i] = PML_C11(r(i), sigma_r_tilde_E[i]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CETHR_TILDE_00[i] = PML_C00(sigma_r_E[i]);
        CETHR_TILDE_01[i] = PML_C01(sigma_r_E[i]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CEPHR_10[i] = PML_C10(r(i), sigma_r_tilde_E[i]);
        CEPHR_11[i] = PML_C11(r(i), sigma_r_tilde_E[i]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CEPHR_TILDE_00[i] = PML_C00(sigma_r_E[i]);
        CEPHR_TILDE_01[i] = PML_C01(sigma_r_E[i]);
    }

    for(int j = 1; j <= Nth - 1; j++){
        CEPHTH_00[j] = PML_C00(sigma_th_E[j]);
        CEPHTH_01[j] = PML_C01(sigma_th_E[j]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CETH_TILDE_00[i] = dt * sigma_r_tilde_E[i] / 2.0;
        CETH_TILDE_01[i] = dt * sigma_r_tilde_E[i] / 2.0;
    }

    for(int i = 1; i <= Nr - 1; i++){
        CEPH_TILDE_00[i] = dt * sigma_r_tilde_E[i] / 2.0;
        CEPH_TILDE_01[i] = dt * sigma_r_tilde_E[i] / 2.0;
    }

    for(int j = 0; j < Nth; j++){
        CHRTH1_00[j] = PML_C00(sigma_th_H[j]);
        CHRTH1_01[j] = PML_C01(sigma_th_H[j]);
    }

    for(int k = 0; k < Nph; k++){
        CHRPH_00[k] = PML_C00(sigma_ph_H[k]);
        CHRPH_01[k] = PML_C01(sigma_ph_H[k]);
    }

    for(int k = 0; k < Nph; k++){
        CHTHPH_00[k] = PML_C00(sigma_ph_H[k]);
        CHTHPH_01[k] = PML_C01(sigma_ph_H[k]);
    }

    for(int i = 0; i < Nr; i++){
        CHTHR_TILDE_00[i] = PML_C00(sigma_r_H[i]);
        CHTHR_TILDE_01[i] = PML_C01(sigma_r_H[i]);
    }

    for(int i = 0; i < Nr; i++){
        CHTHR_10[i] = PML_C10(r(i+0.5), sigma_r_tilde_H[i]);
        CHTHR_11[i] = PML_C11(r(i+0.5), sigma_r_tilde_H[i]);
    }

    for(int j = 0; j < Nth; j++){
        CHPHTH_00[j] = PML_C00(sigma_th_H[j]);
        CHPHTH_01[j] = PML_C01(sigma_th_H[j]);
    }

    for(int i = 0; i < Nr; i++){
        CHPHR_TILDE_00[i] = PML_C00(sigma_r_H[i]);
        CHPHR_TILDE_01[i] = PML_C01(sigma_r_H[i]);
    }

    for(int i = 0; i < Nr; i++){
        CHPHR_10[i] = PML_C10(r(i+0.5), sigma_r_tilde_H[i]);
        CHPHR_11[i] = PML_C11(r(i+0.5), sigma_r_tilde_H[i]);
    }

    for(int i = 0; i < Nr; i++){
        CHTH_TILDE[i] = sigma_r_tilde_H[i] * dt / 2.0;
        CHPH_TILDE[i] = sigma_r_tilde_H[i] * dt / 2.0;
    }

    std::ofstream ofs_C("./data/" + global_dirName + "/C.dat");
    ofs_C << "#CERTH1_00, CERTH1_01, CERPH_00, CERPH_01, CETHPH_00, CETHPH_01, CETHR_10, CETHR_11, CETHR_TILDE_00, CETHR_TILDE_01, CEPHR_10, CEPHR_11, CHTHR_10, CHTHR_11" << std::endl;
    for(int i = 0; i < Nr; i++){
        ofs_C << i << " " << CERTH1_00[i] << " " << CERTH1_01[i] << " " << CERPH_00[i] << " " << CERPH_01[i] << " " 
                << CETHPH_00[i] << " " << CETHPH_01[i] << " " << CETHR_10[i] << " " << CETHR_11[i] << " "
                << CETHR_TILDE_00[i] << " " << CETHR_TILDE_01[i] << " " << CEPHR_10[i] << " " << CEPHR_11[i] 
                << " " << CHTHR_10[i] << " " << CHTHR_11[i]  << std::endl;

    }

    // exit(0);

    delete [] sigma_r_tilde_E;
    delete [] sigma_r_tilde_H;
    delete [] sigma_r_E;
    delete [] sigma_r_H;
    delete [] sigma_th_E;
    delete [] sigma_th_H;
    delete [] sigma_ph_E;
    delete [] sigma_ph_H;

}