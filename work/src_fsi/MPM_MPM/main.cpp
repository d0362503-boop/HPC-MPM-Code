#include <mpi.h>

#include <cmath>
#include <iomanip>
#include <iostream>

#include "module/contact.h"
#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solver/crsmat.h"
#include "work/src_fsi/MPM_MPM/block_fsi.h"

using namespace implicitmpm;
using namespace stabilizedmpm;
using namespace mpmmpmblockfsi;

void mpmmpmblockfsi::MPMBlockFSI() {
    MPMMPMBlockFSI fsi;

    fsi.DataInput();

    istep = ista - 1;
    int iview = istep / iout;
    real_time = dt * double(istep);

    if (rstflag == 1 || rstflag == 3) {
        fsi.solid_.RestartInput();
        fsi.fluid_.RestartInput();
    }

    BuildMesh();

    BuildControlPoint();

    fsi.solid_.SM_.BuildCrsMat(9);

    fsi.fluid_.NS_.BuildCrsMat(16);

    ComputeNodalVol();

    fsi.solid_.MeshPointLinklist();

    fsi.solid_.DetermineRigidBC();

    fsi.solid_.OutputPointDataVTKHDF(iview, istep);

    fsi.fluid_.OutputPointDataVTKHDF(iview, istep);

    for (istep = ista; istep <= iend; istep++) {
        // -----------------------------------------------
        real_time = dt * double(istep);
        // if (real_time < 2.0e0) {
        //     facl = (1.0e0 - std::cos(M_PI * real_time / 2.0e0)) / 2.0e0;
        // } else {
        //     facl = 1.0e0;
        // }
        if (istep <= nlstep) {
            facl = dlstep * double(istep);
        } else {
            facl = 1.0e0;
        }
        // -----------------------------------------------

        fsi.solid_.MeshPointLinklist();

        fsi.fluid_.MeshPointLinklist();

        fsi.solid_.Particle2Node();

        fsi.fluid_.Particle2Node();

        fsi.SolveFSISystem(); // --- Strong coupling: block iteration ---

        fsi.solid_.Node2Particle();

        fsi.fluid_.Node2Particle();

        fsi.solid_.MoveParticle();

        fsi.fluid_.MoveParticle();

        if (istep % iout == 0) {
            iview++;
            fsi.solid_.OutputPointDataVTKHDF(iview, istep);
            fsi.fluid_.OutputPointDataVTKHDF(iview, istep);

            if (rstflag == 2 || rstflag == 3) {
                fsi.solid_.RestartOutput();
                fsi.fluid_.RestartOutput();
            }
        }

        if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, istep); }
    }

    return;
}
