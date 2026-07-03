#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include "../../module/solid/explicit/explicit_mpm_solid.h"
#include "../../module/solid/solid_material_point.h"

using namespace explicitmpm;

void Solid_Explicit_ULMPM() {

    ExplicitSolidMPM sp;

    sp.DataInput();

    if (rstflag == 1 || rstflag == 3) sp.RestartInput();

    BuildMesh();

    BuildControlPoint();

    MakNodalVol();

    istep = ista - 1;
    int iview = istep / iout;
    real_time = dt * double(istep);

    sp.OutputPointDataVTKHDF(iview, istep);

    for (istep = ista; istep <= iend; istep++) {

        // -----------------------------------------------
        real_time = dt * double(istep);
        if (istep <= nlstep) {
            facl = dlstep * double(istep);
        } else {
            facl = 1.0e0;
        }
        // -----------------------------------------------

        sp.MeshPointLinklist();

        sp.Particle2Node();

        sp.SolveSolid();

        sp.Node2Particle();

        if (nprocs != 1) sp.Moveparticle();

        if (istep / iout * iout == istep) {
            iview++;
            sp.OutputPointDataVTKHDF(iview, istep);
        }

        if (istep % 100 == 0 && myrank == 0) {
            std::cout << " ===== Now Computing ===== " << "\n";
            std::cout << "istep:" << std::setw(10) << istep << "\n";
            std::cout << "time:" << std::setw(15) << real_time << std::setw(10) << "iview:" << std::setw(10) << iview << "\n";
            std::cout << " ========================== " << "\n";
        }
    }

    if (rstflag == 2 || rstflag == 3) sp.RestartOutput();

    return;
}
