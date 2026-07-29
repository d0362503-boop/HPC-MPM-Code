#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include "../../module/solid/implicit/implicit_mpm_solid.h"
#include "../../module/solid/solid_material_point.h"
#include "../../module/solver/crsmat.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>

using namespace implicitmpm;

void Solid_implicit_ULMPM() {

    ImplicitSolidMPM sp;

    sp.DataInput();

    if (rstflag == 1 || rstflag == 3) sp.RestartInput();

    BuildMesh();

    BuildControlPoint();

    sp.SM_.BuildCrsMat(9);

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

        sp.MoveParticle();

        if (istep / iout * iout == istep) {
            iview++;
            sp.OutputPointDataVTKHDF(iview, istep);
        }

        if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, istep); }
    }

    if (rstflag == 2 || rstflag == 3) sp.RestartOutput();

    return;
}