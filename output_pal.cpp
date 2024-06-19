#include "fdtd3d.h"
#include <iostream>
#include <fstream>
#include <string>

void output_pal(){
    std::ofstream ofs_pal("./data/" + global_dirName + "/pal.dat",std::ios::app);

    ofs_pal << "真空中　PMLあり" << std::endl;
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

    ofs_pal << "dt = " << dt << std::endl;
    ofs_pal << "Tmax = " << Tmax << std::endl;


    ofs_pal << "source_r = " << int((source_r) / dr) << std::endl;
    ofs_pal << "source_th = " << int((source_th) / Rdth) << std::endl;
    ofs_pal << "source_ph = " << int((source_ph) / Rdph) << std::endl;

    ofs_pal << "PML_L = " << PML_L << std::endl;
    ofs_pal << "PML_R = " << PML_R << std::endl;
    ofs_pal << "PML_M = " << PML_M << std::endl;

    ofs_pal.close();
}