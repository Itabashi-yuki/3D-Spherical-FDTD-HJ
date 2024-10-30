#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <complex>

#include <eigen3/Eigen/Dense>
#include <omp.h>

int main(){
    double ***Er = allocate_3d(Nr, Nth + 1, Nph + 1, 0.0);
    double ****Eth = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);   
    double ****Eph = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    // double ****Erth1 = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    // double ****Erth2 = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    // double ****Erph = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    // double ****Ethr = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    // double ****Ethph = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    // double ****Ephr = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    // double ****Ephth = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);   
    // double ****Ethr_tilde = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    // double ****Ephr_tilde = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Eth_tilde = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);   
    double ****Eph_tilde = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    
    double ****Dr = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Dth = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Dph = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Drth1 = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Drth2 = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Drph = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Dthr = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Dthph = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Dphr = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Dphth = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Dthr_tilde = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Dphr_tilde = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);

    double ***Hr = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ****Hth = allocate_4d(2, Nr, Nth + 1, Nph, 0.0);
    double ****Hph = allocate_4d(2, Nr, Nth+1, Nph + 1, 0.0);
    double ***Hrth1 = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ***Hrth2 = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ***Hrph = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ***Hthr = allocate_3d(Nr, Nth + 1, Nph, 0.0);
    double ***Hthph = allocate_3d(Nr, Nth + 1, Nph, 0.0);
    double ***Hphr = allocate_3d(Nr, Nth, Nph + 1, 0.0);
    double ***Hphth = allocate_3d(Nr, Nth, Nph + 1, 0.0);
    double ****Hthr_tilde = allocate_4d(2, Nr, Nth + 1, Nph, 0.0);
    double ****Hphr_tilde = allocate_4d(2, Nr, Nth, Nph + 1, 0.0);
    double ***Hth_tilde = allocate_3d(Nr, Nth + 1, Nph, 0.0);
    double ***Hph_tilde = allocate_3d(Nr, Nth, Nph + 1, 0.0);
    
    double ****Jr = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Jth = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Jph = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);

    // double *CERTH1_00 = allocate_1d(Nth, 0.0);
    // double *CERTH1_01 = allocate_1d(Nth, 0.0);
    // double *CERPH_00 = allocate_1d(Nph, 0.0);
    // double *CERPH_01 = allocate_1d(Nph, 0.0);
    // double *CETHR_10 = allocate_1d(Nr, 0.0);
    // double *CETHR_11 = allocate_1d(Nr, 0.0);
    // double *CETHPH_00 = allocate_1d(Nph, 0.0);
    // double *CETHPH_01 = allocate_1d(Nph, 0.0);
    // double *CETHR_TILDE_00 = allocate_1d(Nr, 0.0);
    // double *CETHR_TILDE_01 = allocate_1d(Nr, 0.0);
    // double *CEPHR_TILDE_00 = allocate_1d(Nr, 0.0);
    // double *CEPHR_TILDE_01 = allocate_1d(Nr, 0.0);
    // double *CEPHR_10 = allocate_1d(Nr, 0.0);
    // double *CEPHR_11 = allocate_1d(Nr, 0.0);
    // double *CEPHTH_00 = allocate_1d(Nth, 0.0);
    // double *CEPHTH_01 = allocate_1d(Nth, 0.0);
    double *CETH_TILDE_00 = allocate_1d(Nr, 0.0);
    double *CETH_TILDE_01 = allocate_1d(Nr, 0.0);
    double *CEPH_TILDE_00 = allocate_1d(Nr, 0.0);
    double *CEPH_TILDE_01 = allocate_1d(Nr, 0.0);

    double *CDRTH1_00 = allocate_1d(Nth, 0.0);
    double *CDRTH1_01 = allocate_1d(Nth, 0.0);
    double *CDRPH_00 = allocate_1d(Nph, 0.0);
    double *CDRPH_01 = allocate_1d(Nph, 0.0);
    double *CDTHPH_00 = allocate_1d(Nph, 0.0);
    double *CDTHPH_01 = allocate_1d(Nph, 0.0);
    double *CDTHR_10 = allocate_1d(Nr, 0.0);
    double *CDTHR_11 = allocate_1d(Nr, 0.0);
    double *CDTHR_TILDE_00 = allocate_1d(Nr, 0.0);
    double *CDTHR_TILDE_01 = allocate_1d(Nr, 0.0);
    double *CDPHR_TILDE_00 = allocate_1d(Nr, 0.0);
    double *CDPHR_TILDE_01 = allocate_1d(Nr, 0.0);
    double *CDPHR_10 = allocate_1d(Nr, 0.0);
    double *CDPHR_11 = allocate_1d(Nr, 0.0);
    double *CDPHTH_00 = allocate_1d(Nth, 0.0);
    double *CDPHTH_01 = allocate_1d(Nth, 0.0);

    double *CHRTH1_00 = allocate_1d(Nth, 0.0);
    double *CHRTH1_01 = allocate_1d(Nth, 0.0);
    double *CHRPH_00 = allocate_1d(Nph, 0.0);
    double *CHRPH_01 = allocate_1d(Nph, 0.0);
    double *CHTHPH_00 = allocate_1d(Nph, 0.0);
    double *CHTHPH_01 = allocate_1d(Nph, 0.0);
    double *CHTHR_TILDE_00 = allocate_1d(Nr, 0.0);
    double *CHTHR_TILDE_01 = allocate_1d(Nr, 0.0);
    double *CHTHR_10 = allocate_1d(Nr, 0.0);
    double *CHTHR_11 = allocate_1d(Nr, 0.0);
    double *CHPHTH_00 = allocate_1d(Nth, 0.0);
    double *CHPHTH_01= allocate_1d(Nth, 0.0);
    double *CHPHR_TILDE_00 = allocate_1d(Nr, 0.0);
    double *CHPHR_TILDE_01 = allocate_1d(Nr, 0.0);
    double *CHPHR_10 = allocate_1d(Nr, 0.0);
    double *CHPHR_11 = allocate_1d(Nr, 0.0);
    double *CHTH_TILDE = allocate_1d(Nr, 0.0);
    double *CHPH_TILDE = allocate_1d(Nr, 0.0);

    Eigen::Matrix3d ***S = new Eigen::Matrix3d **[Nr+1];
    Eigen::Matrix3d ***B = new Eigen::Matrix3d **[Nr+1];

    for(int i = 0; i < Nr+1; i++){
            S[i] = new Eigen::Matrix3d *[Nth+1];
            B[i] = new Eigen::Matrix3d *[Nth+1];
        for(int j = 0; j < Nth + 1; j++){
            S[i][j] = new Eigen::Matrix3d[Nph+1];
            B[i][j] = new Eigen::Matrix3d[Nph+1];
            for(int k = 0; k < Nph + 1; k++){
                S[i][j][k] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                B[i][j][k] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;        
            }
        }
    }

    double ****check = allocate_4d(2, Nr+1, Nth+1, Nph+1, 0.0);
    make_dir();
    initialize_PML(CDRTH1_00, CDRTH1_01, CDRPH_00, CDRPH_01, CDTHPH_00, CDTHPH_01, CDTHR_10, CDTHR_11, CDTHR_TILDE_00, CDTHR_TILDE_01,
                     CDPHR_10, CDPHR_11, CDPHR_TILDE_00, CDPHR_TILDE_01, CDPHTH_00, CDPHTH_01, CETH_TILDE_00, CETH_TILDE_01, CEPH_TILDE_00,
                      CEPH_TILDE_01, CHRTH1_00, CHRTH1_01,CHRPH_00, CHRPH_01, CHTHPH_00, CHTHPH_01, CHTHR_TILDE_00,CHTHR_TILDE_01, CHTHR_10,
                       CHTHR_11, CHPHTH_00,CHPHTH_01, CHPHR_TILDE_00, CHPHR_TILDE_01, CHPHR_10, CHPHR_11, CHTH_TILDE, CHPH_TILDE);
    initialize_Plasma(S, B);
    // exit(0);
    // std::ofstream ofs("./data/" + global_dirName + "/Coefficient.dat");
    // for(int i = 0; i < Nr; i++){
    //     ofs << CETHR_10[i] << " " << CETHR_11[i] << " " << CETHR_TILDE_00[i] << " " << CETHR_TILDE_01[i] <<  " " << CEPHR_10[i] << " " << CEPHR_11[i] << std::endl;
    // }
    // int n0 = cal_obs_n0();
    double st, en;
    int n0 = 1;
    // std::cout << Nt << std::endl;
    // exit(0);
    std::ofstream ofs_div_time( PATH + "data/" + global_dirName +"div_time_dt_"+ std::to_string(exp_Ne) + ".dat", std::ios::app);
    std::ofstream ofs_obs( PATH + "data/" + global_dirName +"obs.dat",std::ios::app);
    // ofs_obs << "r方向 + 5km, th方向 + 5km, ph方向 + 5km " << std::endl;

    std::complex <double> zj { 0., 1. };
    std::complex <double> *Er0 = new std::complex <double> [Nph+1 - PML_L];
    for(int i = 0; i <= Nph - PML_L; i++){
        Er0[i] = std::complex <double> {0., 0.};
    }

    // double Ave = 0.;
    // std::ofstream ofs_source(PATH + "data/" + global_dirName + "source.dat", std::ios::app);
    for(int n = 1; n < 300; n++){
        int NEW = n % 2;
        if(n % 100 == 0){
            std::cout << n << " / " << Nt << std::endl;
        }
        // if(n % 1000 == 0){
        //     std::cout << n << " / " << Nt << std::endl;
        // }
        double t = dt * ( n - 0.5 );
        // st = omp_get_wtime();
        update_Er(Er, Hth, Hph, Jr, check, n);
        update_Eth(Eth, Hr, Hph, Jth, check, n);
        update_Eph(Eph, Hr, Hth, Jph, check, n);
        // std::cout << en - st << std::endl;
        update_Dr_PML(Drth1, Drth2, Drph, Dr, Hr, Hth, Hph, CDRTH1_00, CDRTH1_01, CDRPH_00, CDRPH_01, check, n);
        update_Dth_PML(Dthph, Dthr, Dthr_tilde, Dth, Hr, Hph_tilde, CDTHPH_00, CDTHPH_01,
                         CDTHR_10, CDTHR_11, CDTHR_TILDE_00, CDTHR_TILDE_01, check, n);
        update_Dph_PML(Dphr, Dphr_tilde, Dphth, Dph, Hr, Hth_tilde, CDPHR_10, CDPHR_11,
                         CDPHR_TILDE_00, CDPHR_TILDE_01, CDPHTH_00, CDPHTH_01, check, n);

        
        update_Er_PML(Er, Dr, Hth, Hph, Jr, check, n);
        update_Eth_PML(Eth, Dth, Hr, Hph, Jth, check, n);
        update_Eph_PML(Eph, Dph, Hr, Hth, Jph, check, n);
        update_Eth_tilde(Eth_tilde, Eth, CETH_TILDE_00,check, n);
        update_Eph_tilde(Eph_tilde, Eph, CEPH_TILDE_00, check,  n);

       
        Er[int((source_r) / dr)][int((source_th) / Rdth)][int((source_ph) / Rdph)] -= dt / EPS0 * source_J(t);

        update_Hr(Hr, Eth, Eph, check, n);
        update_Hth(Hth, Er, Eph, check, n);
        update_Hph(Hph, Er, Eth, check, n);


        update_Hr_PML(Hr, Hrth1, Hrth2, Hrph, Eth, Eph,
                         CHRTH1_00, CHRTH1_01, CHRPH_00, CHRPH_01, check, n);
        update_Hth_PML(Hth, Hthr, Hthph, Hthr_tilde, Er, Eph_tilde,
                         CHTHPH_00, CHTHPH_01, CHTHR_TILDE_00, CHTHR_TILDE_01,
                         CHTHR_10, CHTHR_11, check, n);
        update_Hph_PML(Hph, Hphr, Hphth, Hphr_tilde, Er, Eth_tilde,
                         CHPHTH_00, CHPHTH_01, CHPHR_TILDE_00, CHPHR_TILDE_01,
                          CHPHR_10, CHPHR_11, check, n);
        update_Hth_tilde(Hth_tilde, Hth, CHTH_TILDE, check, n);
        update_Hph_tilde(Hph_tilde, Hph, CHPH_TILDE, check,n);
        update_Jr(Jr, Jth, Jph, Er, Eth, Eph, S, B, n);
        // exit(0);
        update_Jth(Jr, Jth, Jph, Er, Eth, Eph, S, B, n);
        update_Jph(Jr, Jth, Jph, Er, Eth, Eph, S, B, n);
        // en = omp_get_wtime();
        // Ave += en - st;
        // std::ofstream ofs("./data/" + global_dirName + "/E_PML_" + std::to_string(n) + ".dat");
        // for(int k = 0; k <= Nph; k+=2){
        //     ofs << k * R0 * dph * 1.0e-3 << " " << Erth1[NEW][Nr / 2][Nth / 2][k] << " " << Erth2[NEW][Nr / 2][Nth / 2][k] << " " << Erph[NEW][Nr / 2][Nth / 2][k] 
        //             << " " << Ethph[NEW][Nr / 2][Nth / 2][k] << " " << Ethr[NEW][Nr / 2][Nth / 2][k] << " " << Ethr_tilde[NEW][Nr / 2][Nth / 2][k] << std::endl;          
        // }
        output_E(Er, Eth, Eph, Hr, Hth, Hph, Jr, Jth, Jph, Dr, Dth, Dph, n, n0);

        if(std::abs(Er[Nr-5][Nth/2][1]) > 10000){
            ofs_div_time << Er[obs_Nr][obs_Nth][obs_Nph] << " " << n * dt << std::endl;
            output_pal();
            std::cout << "発散しました" << std::endl;
            break;
        }
        ofs_obs << n * dt   << " " << Er[obs_Nr][obs_Nth][obs_Nph]
                            << " " << Er[int((Rr_iono_lower + 5.0e3) / dr)][obs_Nth][obs_Nph]
                            << " " << Er[PML_L][Nth/2][Nph-PML_L - 1] << std::endl;
        
        for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
            Er0[k] += Er[PML_L][Nth/2][k] * std::exp( -1.0 * zj * OMG * t ) * dt;
        }

        // ofs_source << n * dt << " " << source_J(t) << std::endl;

    }
    // std::cout << Ave / 10 << std::endl;

    std::ofstream ofs_Er0( PATH + "data/" + global_dirName +"Er_surface.dat");
    for(int k = PML_L + 1; k <= Nph - PML_L - 1; k++){
        ofs_Er0 << k * R0 * dph * 1e-3 << " " << std::abs( Er0[k] ) << " "
            << std::arg( Er0[k] ) << std::endl;
    }

    // std::ofstream ofs_check( PATH + "data/"+ global_dirName + "check.dat");
    // ofs_check << "# n, i, j, k, check" << std::endl;
    //     for(int n = 0; n < 2; n++){
    //         for(int i = 0; i < Nr; i++){
    //             for(int j = 0; j < Nth; j++){
    //                 for(int k = 0; k < Nph; k++){
    //                 ofs_check << n << " " << i << " " << j << " " << k << " " << check[n][i][j][k] << std::endl;
    //                 }
    //                 ofs_check << std::endl;
    //             }
    //             ofs_check << std::endl;
    //         }
    //         ofs_check << std::endl;
    //     }

    // for(int i = 0; i < Nr; i++){
    //     for(int j = 0; j < Nth; j++){
    //         ofs_check << i << " " << j << " " << check[0][i][j][50] << std::endl;
    //     }
    //     ofs_check << std::endl;
    // }

    output_pal();
    
    free_memory3d(Er, Nr, Nth+1);
    free_memory4d(Eth, 2, Nr + 1, Nth);
    free_memory4d(Eph, 2, Nr + 1, Nth + 1);
    free_memory4d(Dr, 2, Nr, Nth + 1);
    free_memory4d(Dth, 2, Nr + 1, Nth);
    free_memory4d(Dph, 2, Nr + 1, Nth + 1);
    // free_memory4d(Erth1, 2, Nr, Nth+1);
    // free_memory4d(Erth2, 2, Nr, Nth+1);
    // free_memory4d(Erph, 2, Nr, Nth+1);
    // free_memory4d(Ethr, 2, Nr+1, Nth);
    // free_memory4d(Ethph, 2, Nr+1, Nth);
    // free_memory4d(Ephr, 2, Nr+1, Nth+1);
    // free_memory4d(Ephth, 2, Nr+1, Nth+1);
    // free_memory4d(Ethr_tilde, 2, Nr+1, Nth);
    // free_memory4d(Ephr_tilde, 2, Nr+1, Nth+1);
    free_memory4d(Eth_tilde, 2, Nr + 1, Nth);
    free_memory4d(Eph_tilde, 2, Nr + 1, Nth + 1);

    free_memory4d(Drth1, 2, Nr, Nth+1);
    free_memory4d(Drth2, 2, Nr, Nth+1);
    free_memory4d(Drph, 2, Nr, Nth+1);
    free_memory4d(Dthr, 2, Nr+1, Nth);
    free_memory4d(Dthph, 2, Nr+1, Nth);
    free_memory4d(Dphr, 2, Nr+1, Nth+1);
    free_memory4d(Dphth, 2, Nr+1, Nth+1);
    free_memory4d(Dthr_tilde, 2, Nr+1, Nth);
    free_memory4d(Dphr_tilde, 2, Nr+1, Nth+1);

    free_memory3d(Hrth1, Nr+1, Nth);
    free_memory3d(Hrth2, Nr+1, Nth);
    free_memory3d(Hrph, Nr+1, Nth);
    free_memory3d(Hr, Nr + 1, Nth);
    free_memory4d(Hth, 2, Nr, Nth + 1);
    free_memory3d(Hth_tilde, Nr, Nth + 1);
    free_memory3d(Hthr, Nr, Nth + 1);
    free_memory4d(Hthr_tilde, 2, Nr, Nth + 1);
    free_memory3d(Hthph, Nr, Nth + 1);
    free_memory4d(Hph, 2, Nr, Nth);
    free_memory3d(Hph_tilde, Nr, Nth);
    free_memory3d(Hphr, Nr, Nth);
    free_memory4d(Hphr_tilde, 2, Nr, Nth);
    free_memory3d(Hphth, Nr, Nth);

    free_memory4d(Jr, 2, Nr, Nth + 1);
    free_memory4d(Jth, 2, Nr + 1, Nth);
    free_memory4d(Jph, 2, Nr + 1, Nth + 1);

    // free_memory3d(S, Nr+1, Nth+1);
    // free_memory3d(B, Nr+1, Nth+1);

    delete [] CDRTH1_00;
    delete [] CDRTH1_01;
    delete [] CDRPH_00;
    delete [] CDRPH_01;
    delete [] CDTHPH_00;
    delete [] CDTHPH_01;
    delete [] CDTHR_10;
    delete [] CDTHR_11;
    delete [] CDTHR_TILDE_00;
    delete [] CDTHR_TILDE_01;
    delete [] CDPHR_TILDE_00;
    delete [] CDPHR_TILDE_01;
    delete [] CDPHR_10;
    delete [] CDPHR_11;
    delete [] CDPHTH_00;
    delete [] CDPHTH_01;
    delete [] CETH_TILDE_00;
    delete [] CETH_TILDE_01;
    delete [] CEPH_TILDE_00;
    delete [] CEPH_TILDE_01;
    delete [] CHRTH1_00;
    delete [] CHRTH1_01;
    delete [] CHRPH_00;
    delete [] CHRPH_01;
    delete [] CHTHPH_00;
    delete [] CHTHPH_01;
    delete [] CHTHR_TILDE_00;
    delete [] CHTHR_TILDE_01;
    delete [] CHTHR_10;
    delete [] CHTHR_11;
    delete [] CHPHTH_00;
    delete [] CHPHTH_01;
    delete [] CHPHR_TILDE_00;
    delete [] CHPHR_TILDE_01;
    delete [] CHPHR_10;
    delete [] CHPHR_11;
    delete [] CHTH_TILDE;
    delete [] CHPH_TILDE;
}