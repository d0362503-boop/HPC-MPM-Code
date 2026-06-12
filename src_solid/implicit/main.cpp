#include <cmath>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "../../module/dataset.h"
#include "../../module/mesh.h"
#include "../../module/restart.h"
#include "../../module/mpi_data.h"
#include "../../module/solver/crsmat.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/solid/solid_material_point.h"
#include "../../module/solid/implicit/implicit_mpm_solid.h"

using namespace ImplicitMPM;

// -----------------------
void datain_para();
void Particle2NodeS();
void solve_linear_system(CrsMat& solid);
void Node2ParticleS();
// -----------------------

void Solid_implicit_ULMPM () {

    datain_para();

    if (rstflag == 1 || rstflag == 3) sp.RestartInput();

    BuildMesh();

    BuildControlPoint();

    sp.solid.BuildCrsMat(9);

    sp.MakNodeVol();

    istep = ista - 1;
    int iview = istep / iout;
    real_time = dt * double(istep);

    sp.OutputPointData(iview, istep);
    
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

        sp.Particle2Node();

        sp.SolveSolid();

        sp.Node2Particle();

        if (nprocs != 1) sp.Moveparticle();

        if (istep/iout*iout == istep) {
            iview++;
            sp.OutputPointData(iview, istep);
        }

        if (istep%100 == 0 && myrank == 0) {
            std::cout<< " ===== Now Computing ===== " << "\n";
            std::cout<< "istep:" << std::setw(10) << istep << "\n"; 
            std::cout<< "time:" << std::setw(15) << real_time 
                 << std::setw(10) << "iview:" << std::setw(10) << iview << "\n";
            std::cout<< "========================== " << "\n";
        }
    }

    if (rstflag == 2 || rstflag == 3) sp.RestartOutput();

    return;
}