#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "work/src_fsi/MPM_FEM/block_fsi.h"

using namespace mpmfemblockfsi;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        Grid grid;          // rank grid
        RunParams params;   // run parameters
        TimeState ts;       // time-step state

        MPMFEMBlockFSI fsi(grid, params);

        fsi.DataInput(grid, params);

        ts.istep = params.ista - 1;
        int iview = ts.istep / params.iout;
        ts.real_time = params.dt * double(ts.istep);

        if (params.rstflag == 1 || params.rstflag == 3) {
            fsi.solid_.RestartInput();
            fsi.fluid_.RestartInput();
        }

        BuildMesh(grid.mesh);

        BuildControlPoint(grid.mesh, grid.cmesh);

        fsi.fluid_.BuildSubcellPoints();

        fsi.fluid_.MeshPointLinklist();

        fsi.fluid_.Particle2NodePhi();

        double vol_l0 = fsi.fluid_.CalLiquidVol();

        fsi.solid_.SM_.BuildCrsMat(9);

        fsi.fluid_.NS_.BuildCrsMat(16);

        ComputeNodalVol(params, grid);

        fsi.solid_.MeshPointLinklist();

        fsi.solid_.DetermineRigidBC();

        fsi.fluid_.Cp2NodeVTK();

        fsi.solid_.OutputPointDataVTKHDF(iview, ts);

        fsi.fluid_.OutputMeshDataVTKHDF(iview, ts);

        for (ts.istep = params.ista; ts.istep <= params.iend; ts.istep++) {
            // -----------------------------------------------
            ts.real_time = params.dt * double(ts.istep);
            // if (real_time < 2.0e0) {
            //     facl = (1.0e0 - std::cos(M_PI * real_time / 2.0e0)) / 2.0e0;
            // } else {
            //     facl = 1.0e0;
            // }
            if (ts.istep <= params.nlstep) {
                ts.facl = params.dlstep * double(ts.istep);
            } else {
                ts.facl = 1.0e0;
            }
            // -----------------------------------------------

            fsi.fluid_.SetPFDomain();

            fsi.solid_.MeshPointLinklist();

            fsi.solid_.Particle2Node(ts);

            fsi.SolveFSISystem(ts); // --- Strong coupling: block iteration ---

            fsi.solid_.Node2Particle(ts);

            fsi.solid_.MoveParticle();

            if (ts.istep % params.iout == 0) {
                iview++;
                fsi.solid_.OutputPointDataVTKHDF(iview, ts);
                fsi.fluid_.Cp2NodeVTK();
                fsi.fluid_.OutputMeshDataVTKHDF(iview, ts);

                if (params.rstflag == 2 || params.rstflag == 3) {
                    fsi.solid_.RestartOutput();
                    fsi.fluid_.RestartOutput();
                }
            }

            if (ts.istep % 100 == 0 && myrank == 0) { OutputMessage(iview, ts); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
