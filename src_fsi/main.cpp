#include <mpi.h>

#include <cmath>
#include <iomanip>
#include <iostream>

#include "../module/contact.h"
#include "../module/data_io.h"
#include "../module/dataset.h"
#include "../module/fluid/FEM/stabilized_fem.h"
#include "../module/material_point.h"
#include "../module/mesh.h"
#include "../module/mpi_data.h"
#include "../module/solid/implicit/implicit_mpm_solid.h"
#include "../module/solver/crsmat.h"
#include "block_fsi.h"

using namespace implicitmpm;
using namespace stabilizedfem;

void MPMBlockFSI() {
    BlockFSI fsi;

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

    fsi.fluid_.BuildGaussianPoint();

    fsi.fluid_.MeshPointLinklist();

    fsi.fluid_.Particle2NodePhi();

    double vol_l0 = fsi.fluid_.CalLiquidVol();

    fsi.solid_.SM_.BuildCrsMat(9);

    fsi.fluid_.NS_.BuildCrsMat(16);

    MakNodalVol();

    fsi.solid_.MeshPointLinklist();

    fsi.solid_.DetermineRigidBC();

    fsi.fluid_.Cp2NodeVTK();

    fsi.solid_.OutputPointData(iview, istep);

    fsi.fluid_.OutputMeshData(iview, istep);

    for (istep = ista; istep <= iend; istep++) {
        // -----------------------------------------------
        real_time = dt * double(istep);
        if (real_time < 2.0e0) {
            facl = (1.0e0 - cos(M_PI * real_time / 2.0e0)) / 2.0e0;
        } else {
            facl = 1.0e0;
        }
        // if (istep <= nlstep) {
        //     facl = dlstep * double(istep);
        // } else {
        //     facl = 1.0e0;
        // }
        // -----------------------------------------------

        fsi.fluid_.SetPFDomain();

        fsi.solid_.MeshPointLinklist();

        fsi.solid_.Particle2Node();

        fsi.SolveFSISystem(); // --- Strong coupling: block iteration ---

        fsi.solid_.Node2Particle();

        if (nprocs != 1) fsi.solid_.Moveparticle();

        if (istep % iout == 0) {
            iview++;
            fsi.solid_.OutputPointData(iview, istep);
            fsi.fluid_.Cp2NodeVTK();
            fsi.fluid_.OutputMeshData(iview, istep);

            if (rstflag == 2 || rstflag == 3) {
                fsi.solid_.RestartOutput();
                fsi.fluid_.RestartOutput();
            }
        }

        if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, istep); }
    }

    return;
}
