#include <mpi.h>

#include <cmath>
#include <iomanip>
#include <iostream>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solver/crsmat.h"
#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace implicitmpm;
using namespace stabilizedmpm;
using namespace mpm_mpm_monolithic_fsi;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        MPMMPMMonolithicFSI fsi;

        fsi.DataInput();

        istep = ista - 1;
        int iview = istep / iout;
        real_time = dt * double(istep);

        if (rstflag == 1 || rstflag == 3) {
            fsi.fluid_.RestartInput();
            fsi.solid_.RestartInput();
        }

        BuildMesh();

        BuildControlPoint();

        fsi.fsi_sys.BuildCrsMat(MPMMPMMonolithicFSI::kStoredBlocks);

        ComputeNodalVol();

        fsi.solid_.MeshPointLinklist();

        fsi.solid_.DetermineRigidBC();

        fsi.solid_.OutputPointDataVTKHDF(iview, istep);

        fsi.fluid_.OutputPointDataVTKHDF(iview, istep);

        for (istep = ista; istep <= iend; istep++) {
            // -----------------------------------------------
            real_time = dt * double(istep);
            double s_t;
            if (real_time < 0.51e0) {
                s_t = std::pow(real_time, 2) / (2.0e0 * 0.51e0);
            } else {
                s_t = real_time - 0.51e0 / 2.0e0;
            }
            double theta = M_PI / 45.0e0 * std::sin(2.0e0 * M_PI / 1.211e0 * s_t);
            bb[0] = -9.81e0 * std::sin(theta);
            bb[2] = -9.81e0 * std::cos(theta);

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

            fsi.SolveFSISystem(); // --- Monolithic Newton iteration ---

            fsi.solid_.Node2Particle();

            fsi.fluid_.Node2Particle();

            fsi.solid_.MoveParticle();

            fsi.CorrectFluidPenetration();

            fsi.fluid_.MoveParticle();

            if (istep % iout == 0) {
                iview++;
                fsi.solid_.OutputPointDataVTKHDF(iview, istep);
                fsi.fluid_.OutputPointDataVTKHDF(iview, istep);

                if (rstflag == 2 || rstflag == 3) {
                    fsi.fluid_.RestartOutput();
                    fsi.solid_.RestartOutput();
                }
            }

            if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, istep); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
