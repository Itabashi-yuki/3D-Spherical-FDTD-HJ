#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <string>

void output_pal(){
    std::ofstream ofs_pal( PATH + "data/" + global_dirName + "pal.dat",std::ios::app);

    ofs_pal << "プラズマ領域あり　PMLは真空中" << std::endl;
    ofs_pal << std::endl;

    ofs_pal << "Rr = " << Rr << std::endl;
    ofs_pal << "Rth = " << Rth << std::endl;
    ofs_pal << "Rph = " << Rph << std::endl;
    ofs_pal << "Rdth = " << Rdth << std::endl;
    ofs_pal << "Rdph = " << Rdph << std::endl;
    ofs_pal << "dr = " << dr << std::endl;
    ofs_pal << "dth = " << dth << std::endl;
    ofs_pal << "dph = " << dph << std::endl;
    ofs_pal << "thR = " << thR << std::endl;

    ofs_pal << "xi = " << xi << std::endl;
    ofs_pal << "dt = " << dt << std::endl;
    ofs_pal << "Tmax = " << Tmax << std::endl;


    ofs_pal << "source_r = " << int((source_r) / dr) << std::endl;
    ofs_pal << "source_th = " << int((source_th) / Rdth) << std::endl;
    ofs_pal << "source_ph = " << int((source_ph) / Rdph) << std::endl;

    ofs_pal << "obs_r = " << obs_r << std::endl;
    ofs_pal << "obs_th = " << obs_th << std::endl;
    ofs_pal << "obs_ph = " << obs_ph << std::endl;

    ofs_pal << "-------PMLパラメタ-------" << std::endl;
    ofs_pal << "PML_L = " << PML_L << std::endl;
    ofs_pal << "PML_R = " << PML_R << std::endl;
    ofs_pal << "PML_M = " << PML_M << std::endl;

    ofs_pal << "-------プラズマパラメタ-------" << std::endl;
    ofs_pal << "Rr_iono_lower = " << Rr_iono_lower << std::endl;
    ofs_pal << "Rr_iono_upper = " << Rr_iono_upper << std::endl;
    ofs_pal << "Rth_iono_lower = " << Rth_iono_lower << std::endl;
    ofs_pal << "Rth_iono_upper = " << Rth_iono_upper << std::endl;
    ofs_pal << "Rph_iono_lower = " << Rph_iono_lower << std::endl;
    ofs_pal << "Rph_iono_upper = " << Rph_iono_upper << std::endl;
    ofs_pal << "B0 = " << B0 << std::endl;
    ofs_pal << "Ne = " << cal_Ne(exp_Ne) << std::endl;
    ofs_pal << "nu = " << cal_nu(exp_nu) << std::endl;
    ofs_pal << "THETA = " << THETA << std::endl; 
    ofs_pal << "PHI = " << PHI << std::endl; 

    ofs_pal << "-------大地インピーダンスパラメタ-------" << std::endl;
    ofs_pal << "sigma = " << sig_s(0,0) << std::endl;

    ofs_pal.close();
}