#include <cmath>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        Grid grid;          // rank grid
        RunParams params;   // run parameters
        TimeState ts;       // time-step state

        MPMMPMMonolithicFSI fsi(grid, params);

        fsi.DataInput(grid, params);

        ts.istep = params.ista - 1;
        int iview = ts.istep / params.iout;
        ts.real_time = params.dt * double(ts.istep);

        if (params.rstflag == 1 || params.rstflag == 3) {
            fsi.fluid_.RestartInput();
            fsi.solid_.RestartInput();
        }

        BuildMesh(grid.mesh);

        BuildControlPoint(grid.mesh, grid.cmesh);

        fsi.fsi_sys.BuildCrsMat(MPMMPMMonolithicFSI::kStoredBlocks);

        ComputeNodalVol(params, grid);

        fsi.solid_.MeshPointLinklist();

        fsi.solid_.DetermineRigidBC();

        fsi.solid_.OutputPointDataVTKHDF(iview, ts);

        fsi.fluid_.OutputPointDataVTKHDF(iview, ts);

        for (ts.istep = params.ista; ts.istep <= params.iend; ts.istep++) {
            // -----------------------------------------------
            ts.real_time = params.dt * double(ts.istep);
            // double s_t;
            // if (real_time < 0.51e0) {
            //     s_t = std::pow(real_time, 2) / (2.0e0 * 0.51e0);
            // } else {
            //     s_t = real_time - 0.51e0 / 2.0e0;
            // }
            // double theta = M_PI / 45.0e0 * std::sin(2.0e0 * M_PI / 1.211e0 * s_t);
            // bb[0] = -9.81e0 * std::sin(theta);
            // bb[2] = -9.81e0 * std::cos(theta);

            if (ts.istep <= params.nlstep) {
                ts.facl = params.dlstep * double(ts.istep);
            } else {
                ts.facl = 1.0e0;
            }
            // -----------------------------------------------

            fsi.solid_.MeshPointLinklist();

            fsi.fluid_.MeshPointLinklist();

            fsi.solid_.Particle2Node(ts);

            fsi.fluid_.Particle2Node(ts);

            fsi.SolveFSISystem(ts); // --- Monolithic Newton iteration ---

            fsi.solid_.Node2Particle(ts);

            fsi.fluid_.Node2Particle(ts);

            fsi.solid_.MoveParticle();

            fsi.CorrectFluidPenetration();

            fsi.fluid_.MoveParticle();

            if (ts.istep % params.iout == 0) {
                iview++;
                fsi.solid_.OutputPointDataVTKHDF(iview, ts);
                fsi.fluid_.OutputPointDataVTKHDF(iview, ts);

                if (params.rstflag == 2 || params.rstflag == 3) {
                    fsi.fluid_.RestartOutput();
                    fsi.solid_.RestartOutput();
                }
            }

            if (ts.istep % 100 == 0 && myrank == 0) { OutputMessage(iview, ts); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
