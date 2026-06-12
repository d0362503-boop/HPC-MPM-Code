#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/fluid/MPM/stabilized_mpm.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include "../../module/solver/crsmat.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>

using namespace stabilizedmpm;

void StabilizedMixedMPM() {

    StabilizedMPM wp;

    wp.DataInput();

    if (rstflag == 1 || rstflag == 3) wp.RestartInput();

    BuildMesh();

    BuildControlPoint();

    wp.NS_.BuildCrsMat(16);

    MakNodalVol();

    istep = ista - 1;
    int iview = istep / iout;
    real_time = dt * double(istep);

    wp.OutputPointDataVTKHDF(iview, istep);

    for (istep = ista; istep <= iend; istep++) {

        // -----------------------------------------------
        real_time = dt * double(istep);
        if (istep <= nlstep) {
            facl = dlstep * double(istep);
        } else {
            facl = 1.0e0;
        }
        // -----------------------------------------------

        wp.MeshPointLinklist();

        wp.Particle2Node();

        wp.SolveNS();

        wp.Node2Particle();

        if (nprocs != 1) wp.Moveparticle();

        if (istep % iout == 0) {
            iview++;
            wp.OutputPointDataVTKHDF(iview, istep);

            if (rstflag == 2 || rstflag == 3) wp.RestartOutput();
        }

        if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, istep); }
    }

    return;
}