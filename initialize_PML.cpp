#include "fdtd3d.h"
#include <iostream>
#include <fstream>

double PML_sigma_r(double i){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / dr * std::log(PML_R);
    //sigma_max = 0.0;
    return sigma_max * std::pow((PML_L - i) * dr / PML_L / dr, PML_M);
}

double PML_sigma_th(double j){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / (R0 * dth) * std::log(PML_R);
    //sigma_max = 0.0;
    return sigma_max * std::pow((PML_L - j) / PML_L, PML_M);
}

double PML_sigma_ph(double k){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / (R0 * dph) * std::log(PML_R);
    //sigma_max = 0.0;
    return sigma_max * std::pow((PML_L - k) / PML_L, PML_M);
}


double PML_sigma_r_tilde(double i){
    double sigma_max = - (PML_M + 1.0) * C0 / 2.0 / PML_L / dr * std::log(PML_R);
    //sigma_max = 0.0;

    // int N = 10000;
    // double Rmax = PML_L * dr;
    // double Rmin = i * dr;
    // double Dr = (Rmax - Rmin) / N;
    // double sigma_r_tilde = 0.0;
    // for(int n = 0; n <= N; n++){
    //     // double r = Rmin + i * Dr;
    //     sigma_r_tilde += PML_sigma_r(i) * Dr;
    // }

    // return sigma_r_tilde;
    return sigma_max * ( PML_L * dr / (PML_M + 1.0)) * std::pow(( (PML_L - i) / PML_L ), (PML_M + 1.0));
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

void initialize_PML(double *CDRTH1_00, double *CDRTH1_01, double *CDRPH_00, double *CDRPH_01, double *CDTHPH_00, double *CDTHPH_01, double *CDTHR_10, double *CDTHR_11,
                     double *CDTHR_TILDE_00, double *CDTHR_TILDE_01, double *CDPHR_10, double *CDPHR_11, double *CDPHR_TILDE_00, double *CDPHR_TILDE_01, double *CDPHTH_00,
                     double *CDPHTH_01, double *CETH_TILDE_00, double *CETH_TILDE_01,
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

//    std::ofstream ofs_sigma_i( PATH + "data/" + global_dirName + "sigma_i.dat");
//     ofs_sigma_i << "#sigma_r_E, sigma_r_tilde_E, sigma_r_H" << std::endl;
//     for(int i = 0; i < Nr; i++){
//         ofs_sigma_i << i << " " << sigma_r_E[i] << " " << sigma_r_tilde_E[i] << " " << sigma_r_H[i] << " " << sigma_r_tilde_H[i] << std::endl;
//     } 

//     std::ofstream ofs_sigma_j(PATH + "data/" + global_dirName + "sigma_j.dat");
//     ofs_sigma_j << "#sigma_th_E, sigma_th_H" << std::endl;
//     for(int j = 0; j < Nth; j++){
//         ofs_sigma_j << j << " " << sigma_th_E[j] << " " << sigma_th_H[j] << std::endl;
//     }

//     std::ofstream ofs_sigma_k(PATH + "data/" + global_dirName + "sigma_k.dat");
//     ofs_sigma_k << "#sigma_ph_E, sigma_ph_H" << std::endl;
//     for(int k = 0; k < Nph; k++){
//         ofs_sigma_k << k << " " << sigma_ph_E[k] << " " << sigma_ph_H[k] << std::endl;
//     }

    // for(int j = 0; j < Nth; j++){
    //     ofs_sigma << j << " " << sigma_r_E[j] << " " << sigma_th_E[j] << " " << sigma_ph_E[j] << " " << sigma_r_H[j] << " " << sigma_th_H[j] << " " 
    //                 << sigma_ph_H[j] << " " << sigma_r_tilde_E[j] << " " << sigma_r_tilde_H[j] << std::endl;
    // }

    // ofs_sigma_i.close();
    // ofs_sigma_j.close();
    // ofs_sigma_k.close();
    // exit(0);
    for(int j = 1; j <= Nth - 1; j++){
        CDRTH1_00[j] = PML_C00(sigma_th_E[j]);
        CDRTH1_01[j] = PML_C01(sigma_th_E[j]);
    }

    for(int k = 1; k <= Nph - 1; k++){
        CDRPH_00[k] = PML_C00(sigma_ph_E[k]);
        CDRPH_01[k] = PML_C01(sigma_ph_E[k]);
    }

    for(int k = 1; k <= Nph - 1; k++){
        CDTHPH_00[k] = PML_C00(sigma_ph_E[k]);
        CDTHPH_01[k] = PML_C01(sigma_ph_E[k]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CDTHR_10[i] = PML_C10(r(i), sigma_r_tilde_E[i]);
        CDTHR_11[i] = PML_C11(r(i), sigma_r_tilde_E[i]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CDTHR_TILDE_00[i] = PML_C00(sigma_r_E[i]);
        CDTHR_TILDE_01[i] = PML_C01(sigma_r_E[i]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CDPHR_10[i] = PML_C10(r(i), sigma_r_tilde_E[i]);
        CDPHR_11[i] = PML_C11(r(i), sigma_r_tilde_E[i]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CDPHR_TILDE_00[i] = PML_C00(sigma_r_E[i]);
        CDPHR_TILDE_01[i] = PML_C01(sigma_r_E[i]);
    }

    for(int j = 1; j <= Nth - 1; j++){
        CDPHTH_00[j] = PML_C00(sigma_th_E[j]);
        CDPHTH_01[j] = PML_C01(sigma_th_E[j]);
    }

    for(int i = 1; i <= Nr - 1; i++){
        CETH_TILDE_00[i] = dt * sigma_r_tilde_E[i] / 2.0;
    }

    for(int i = 1; i <= Nr - 1; i++){
        CEPH_TILDE_00[i] = dt * sigma_r_tilde_E[i] / 2.0;
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

    // std::ofstream ofs_C_i( PATH + "data/" + global_dirName + "C_i.dat");
    // ofs_C_i << "#CETHR_10, CETHR_11, CETHR_TILDE_00, CETHR_TILDE_01, CEPHR_10, CEPHR_11, CEPHR_TILDE_00, CEPHR_TILDE_01, CETH_TILDE_00, CEPH_TILDE_00, CHTHR_10, CHTHR_11, CHTHR_TILDE_00, CHTHR_TILDE_01, CHPHR_10, CHPHR_11, CHPHR_TILDE_00, CHPHR_TILDE_01, CHTH_TILDE, CHPH_TILDE" << std::endl;
    // for(int i = 0; i < Nr; i++){
    //     ofs_C_i << i << " " << CDTHR_10[i] << " " << CDTHR_11[i] << " " << CDTHR_TILDE_00[i] << " " << CDTHR_TILDE_01[i] << " " << CDPHR_10[i] << " " << CDPHR_11[i] 
    //             << " " << CDPHR_TILDE_00[i] << " " << CDPHR_TILDE_01[i] << " " << CDTHR_TILDE_00[i] << " " << CDPHR_TILDE_00[i] <<  CHTHR_10[i] << " " << CHTHR_11[i] 
    //             << " " << CHTHR_TILDE_00[i] << " " << CHTHR_TILDE_01[i] << " " << CHPHR_10[i] << " " << CHPHR_11[i] << " " << CHPHR_TILDE_00[i] << " " << CHPHR_TILDE_01[i] 
    //             << " " << CHTH_TILDE[i] << " " << CHPH_TILDE[i] << std::endl;
    // }

    // std::ofstream ofs_C_j( PATH + "data/" + global_dirName + "C_j.dat");
    // ofs_C_j << "#CERTH1_00, CERTH1_01, CEPHTH_00, CEPHTH_01, CHRTH1_00, CHRTH1_01, CHPHTH_00, CHPHTH_01" << std::endl;
    // for(int j = 0; j < Nth; j++){
    //     ofs_C_j << j << " " << CDRTH1_00[j] << " " << CDRTH1_01[j] << " " << CDPHTH_00[j] << " " << CDPHTH_01[j] << " " << CHRTH1_00[j] << " " << CHRTH1_01[j]
    //                 << " " << CHPHTH_00[j] << " " << CHPHTH_01[j] << std::endl;
    // } 

    // std::ofstream ofs_C_k( PATH + "data/" + global_dirName + "C_k.dat");
    // ofs_C_k << "#CERPH_00, CERPH_01, CETHPH_00, CETHPH_01, CHRPH_00, CHRPH_01, CHTHPH_00, CHTHPH_01" << std::endl;
    // for(int k = 0; k < Nph; k++){
    //     ofs_C_k << k << " " << CDRPH_00[k] << " " << CDRPH_01[k] << " " << CDTHPH_00[k] << " " << CDTHPH_01[k] << " " << CHRPH_00[k] << " " << CHRPH_01[k]
    //                 << " " << CHTHPH_00[k] << " " << CHTHPH_01[k] << std::endl;
    // }


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
