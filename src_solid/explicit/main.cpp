#include <cmath>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "../../module/dataset.h"
#include "../../module/mesh.h"
#include "../../module/restart.h"
#include "../../module/mpi_data.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"


// -----------------------
void datain_para();
void Particle2NodeS();
void solve_linear_system();
void Node2ParticleS();
// -----------------------

void Solid_Explicit_ULMPM () {

    datain_para();

    if (rstflag == 1 || rstflag == 3) RestartInSp();

    BuildMesh ();

    BuildControlPoint();

    sp.MakNodeVol();

    istep = ista - 1;
    int iview = istep / iout;
    real_time = dt * double(istep);

    output_sp(iview, istep);

    for (istep = ista; istep <= iend; istep++) {

    // -----------------------------------------------
        real_time = dt * double(istep);
        if (istep <= nlstep) {
            facl = dlstep * double(istep);
        }
        else {
            facl = 1.0e0;
        }
    // -----------------------------------------------

        sp.MeshPointLinklist();

        Particle2NodeS();

        solve_linear_system();

        Node2ParticleS();

        if (nprocs != 1) moveparticle_s();

        if (istep/iout*iout == istep) {
            iview++;
            output_sp(iview, istep);
        }

        if (istep%100 == 0 && myrank == 0) {
            std::cout<< " ===== Now Computing ===== " << "\n";
            std::cout<< "istep:" << std::setw(10) << istep << "\n"; 
            std::cout<< "time:" << std::setw(15) << real_time 
                 << std::setw(10) << "iview:" << std::setw(10) << iview << "\n";
            std::cout<< " ========================== " << "\n";
        }
    }

    if (rstflag == 2 || rstflag == 3) RestartOutSp();

    return;
}