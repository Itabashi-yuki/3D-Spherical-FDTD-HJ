#include "fdtd3d.h"
#include <iostream>
#include <fstream>

int main(){
    double ****Erth1 = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Erth2 = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Erph = allocate_4d(2, Nr, Nth + 1, Nph + 1, 0.0);
    double ****Ethph = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Ethr = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Ethr_tilde = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);
    double ****Ephr = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Ephr_tilde = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Ephth = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);   
    double ***Er = allocate_3d(Nr, Nth + 1, Nph + 1, 0.0);
    double ****Eth_tilde = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);   
    double ****Eth = allocate_4d(2, Nr + 1, Nth, Nph + 1, 0.0);   
    double ****Eph_tilde = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ****Eph = allocate_4d(2, Nr + 1, Nth + 1, Nph, 0.0);
    double ***Hrth1 = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ***Hrth2 = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ***Hrph = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ***Hr = allocate_3d(Nr + 1, Nth, Nph, 0.0);
    double ****Hth = allocate_4d(2, Nr, Nth + 1, Nph, 0.0);
    double ***Hth_tilde = allocate_3d(Nr, Nth + 1, Nph, 0.0);
    double ***Hthr = allocate_3d(Nr, Nth + 1, Nph, 0.0);
    double ****Hthr_tilde = allocate_4d(2, Nr, Nth + 1, Nph, 0.0);
    double ***Hthph = allocate_3d(Nr, Nth + 1, Nph, 0.0);
    double ****Hph = allocate_4d(2, Nr, Nth+1, Nph + 1, 0.0);
    double ***Hph_tilde = allocate_3d(Nr, Nth, Nph + 1, 0.0);
    double ***Hphr = allocate_3d(Nr, Nth, Nph + 1, 0.0);
    double ****Hphr_tilde = allocate_4d(2, Nr, Nth, Nph + 1, 0.0);
    double ***Hphth = allocate_3d(Nr, Nth, Nph + 1, 0.0);

    double *CERTH1_00 = allocate_1d(Nr, 0.0);
    double *CERTH1_01 = allocate_1d(Nr, 0.0);
    double *CERPH_00 = allocate_1d(Nr, 0.0);
    double *CERPH_01 = allocate_1d(Nr, 0.0);
    double *CETHPH_00 = allocate_1d(Nth, 0.0);
    double *CETHPH_01 = allocate_1d(Nth, 0.0);
    double *CETHR_10 = allocate_1d(Nth, 0.0);
    double *CETHR_11 = allocate_1d(Nth, 0.0);
    double *CETHR_TILDE_00 = allocate_1d(Nth, 0.0);
    double *CETHR_TILDE_01 = allocate_1d(Nth, 0.0);
    double *CEPHR_TILDE_00 = allocate_1d(Nph, 0.0);
    double *CEPHR_TILDE_01 = allocate_1d(Nph, 0.0);
    double *CEPHR_10 = allocate_1d(Nph, 0.0);
    double *CEPHR_11 = allocate_1d(Nph, 0.0);
    double *CEPHTH_00 = allocate_1d(Nph, 0.0);
    double *CEPHTH_01 = allocate_1d(Nph, 0.0);
    double *CETH_TILDE_00 = allocate_1d(Nth, 0.0);
    double *CETH_TILDE_01 = allocate_1d(Nth, 0.0);
    double *CEPH_TILDE_00 = allocate_1d(Nph, 0.0);
    double *CEPH_TILDE_01 = allocate_1d(Nph, 0.0);
    double *CHRTH1_00 = allocate_1d(Nr, 0.0);
    double *CHRTH1_01 = allocate_1d(Nr, 0.0);
    double *CHRPH_00 = allocate_1d(Nr, 0.0);
    double *CHRPH_01 = allocate_1d(Nr, 0.0);
    double *CHTHPH_00 = allocate_1d(Nth, 0.0);
    double *CHTHPH_01 = allocate_1d(Nth, 0.0);
    double *CHTHR_TILDE_00 = allocate_1d(Nth, 0.0);
    double *CHTHR_TILDE_01 = allocate_1d(Nth, 0.0);
    double *CHTHR_10 = allocate_1d(Nth, 0.0);
    double *CHTHR_11 = allocate_1d(Nth, 0.0);
    double *CHPHTH_00 = allocate_1d(Nph, 0.0);
    double *CHPHTH_01= allocate_1d(Nph, 0.0);
    double *CHPHR_TILDE_00 = allocate_1d(Nph, 0.0);
    double *CHPHR_TILDE_01 = allocate_1d(Nph, 0.0);
    double *CHPHR_10 = allocate_1d(Nph, 0.0);
    double *CHPHR_11 = allocate_1d(Nph, 0.0);
    double *CHTH_TILDE = allocate_1d(Nth, 0.0);
    double *CHPH_TILDE = allocate_1d(Nph, 0.0);

    double ****check = allocate_4d(2, Nr+1, Nth+1, Nph+1, 0.0);

    make_dir();
    initialize_PML(CERTH1_00, CERTH1_01, CERPH_00, CERPH_01, CETHPH_00, CETHPH_01, CETHR_10, CETHR_11, CETHR_TILDE_00, CETHR_TILDE_01,
                     CEPHR_10, CEPHR_11, CEPHR_TILDE_00, CEPHR_TILDE_01, CEPHTH_00, CEPHTH_01, CETH_TILDE_00, CETH_TILDE_01, CEPH_TILDE_00,
                      CEPH_TILDE_01, CHRTH1_00, CHRTH1_01,CHRPH_00, CHRPH_01, CHTHPH_00, CHTHPH_01, CHTHR_TILDE_00,CHPHR_TILDE_01, CHTHR_10,
                       CHTHR_11, CHPHTH_00,CHPHTH_01, CHPHR_TILDE_00, CHPHR_TILDE_01, CHPHR_10, CHPHR_11, CHTH_TILDE, CHPH_TILDE);

    // exit(0);
    // std::ofstream ofs("./data/" + global_dirName + "/Coefficient.dat");
    // for(int i = 0; i < Nr; i++){
    //     ofs << CETHR_10[i] << " " << CETHR_11[i] << " " << CETHR_TILDE_00[i] << " " << CETHR_TILDE_01[i] <<  " " << CEPHR_10[i] << " " << CEPHR_11[i] << std::endl;
    // }
    // exit(0);
    // int n0 = cal_obs_n0();
    int n0 = 2;

    std::ofstream ofs_obs("./data/" + global_dirName +"/obs.dat",std::ios::app);
    for(int n = 1; n < Nt; n++){
        if(n % 10 == 0)
            std::cout << n << " / " << Nt << std::endl;

        double t = dt * ( n - 0.5 );
        update_Er(Er, Hth, Hph, check, n);
        update_Eth(Eth, Hr, Hph, check, n);
        update_Eph(Eph, Hr, Hth, check, n);
        update_Er_PML(Erth1, Erth2, Erph, Er, Hr, Hth, Hph, CERTH1_00, CERTH1_01, CERPH_00, CERPH_01, check, n);
    // exit(0);
        update_Eth_PML(Ethph, Ethr, Ethr_tilde, Eth, Hr, Hph_tilde, CETHPH_00, CETHPH_01,
                         CETHR_10, CETHR_11, CETHR_TILDE_00, CETHR_TILDE_01, check, n);
        update_Eph_PML(Ephr, Ephr_tilde, Ephth, Eph, Hr, Hth_tilde, CEPHR_10, CEPHR_11,
                         CEPHR_TILDE_00, CEPHR_TILDE_01, CEPHTH_00, CEPHTH_01, check, n);
        update_Eth_tilde(Eth_tilde, Eth, CETH_TILDE_00,check, n);
        update_Eph_tilde(Eph_tilde, Eph, CEPH_TILDE_00, check,  n);

        Er[int((source_r) / dr)][int((source_th) / Rdth)][int((source_ph) / Rdph)] -= dt / EPS0 * Jr(t);
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
        update_Hth_tilde(Hth_tilde, Hth, CHTH_TILDE, n);
        update_Hph_tilde(Hph_tilde, Hph, CHPH_TILDE, n);

        // output_E(Er, Eth, Eph, Hr, Hth, Hph, n, n0);

        ofs_obs << n * dt << " " << Er[int(obs_r / dr)][int(obs_th / Rdth)][int(obs_ph / Rdph)] << std::endl;
    }

    std::ofstream ofs_check("./data/"+ global_dirName + "/check.dat");
    ofs_check << "# n, i, j, k, check" << std::endl;
        for(int n = 0; n < 2; n++){
            for(int i = 0; i < Nr; i++){
                for(int j = 0; j < Nth; j++){
                    for(int k = 0; k < Nph; k++){
                    ofs_check << n << " " << i << " " << j << " " << k << " " << check[n][i][j][k] << std::endl;
                    }
                    ofs_check << std::endl;
                }
                ofs_check << std::endl;
            }
            ofs_check << std::endl;
        }

    // for(int i = 0; i < Nr; i++){
    //     for(int j = 0; j < Nth; j++){
    //         ofs_check << i << " " << j << " " << check[0][i][j][50] << std::endl;
    //     }
    //     ofs_check << std::endl;
    // }

    output_pal();
    
    free_memory4d(Erth1, 2, Nr, Nth+1);
    free_memory4d(Erth2, 2, Nr, Nth+1);
    free_memory4d(Erph, 2, Nr, Nth+1);
    free_memory4d(Ethph, 2, Nr+1, Nth);
    free_memory4d(Ethr, 2, Nr+1, Nth);
    free_memory4d(Ethr_tilde, 2, Nr+1, Nth);
    free_memory4d(Ephr, 2, Nr+1, Nth+1);
    free_memory4d(Ephr_tilde, 2, Nr+1, Nth+1);
    free_memory4d(Ephth, 2, Nr+1, Nth+1);
    free_memory3d(Er, Nr, Nth+1);
    free_memory4d(Eth, 2, Nr + 1, Nth);
    free_memory4d(Eph, 2, Nr + 1, Nth + 1);
    free_memory4d(Eth_tilde, 2, Nr + 1, Nth);
    free_memory4d(Eph_tilde, 2, Nr + 1, Nth + 1);

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

    delete [] CERTH1_00;
    delete [] CERTH1_01;
    delete [] CERPH_00;
    delete [] CERPH_01;
    delete [] CETHPH_00;
    delete [] CETHPH_01;
    delete [] CETHR_10;
    delete [] CETHR_11;
    delete [] CETHR_TILDE_00;
    delete [] CETHR_TILDE_01;
    delete [] CEPHR_TILDE_00;
    delete [] CEPHR_TILDE_01;
    delete [] CEPHR_10;
    delete [] CEPHR_11;
    delete [] CEPHTH_00;
    delete [] CEPHTH_01;
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