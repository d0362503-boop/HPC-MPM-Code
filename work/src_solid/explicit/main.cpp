#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/explicit/explicit_mpm_solid.h"

using namespace explicitmpm;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        Grid grid;          // rank grid
        RunParams params;   // run parameters
        TimeState ts;       // time-step state

        ExplicitSolidMPM sp(grid, params);

        sp.DataInput(grid, params);

        if (params.rstflag == 1 || params.rstflag == 3) sp.RestartInput();

        BuildMesh(grid.mesh);

        BuildControlPoint(grid.mesh, grid.cmesh);

        ComputeNodalVol(params, grid);

        ts.istep = params.ista - 1;
        int iview = ts.istep / params.iout;
        ts.real_time = params.dt * double(ts.istep);

        sp.OutputPointDataVTKHDF(iview, ts);

        for (ts.istep = params.ista; ts.istep <= params.iend; ts.istep++) {

            // -----------------------------------------------
            ts.real_time = params.dt * double(ts.istep);
            if (ts.istep <= params.nlstep) {
                ts.facl = params.dlstep * double(ts.istep);
            } else {
                ts.facl = 1.0e0;
            }
            // -----------------------------------------------

            sp.MeshPointLinklist();

            sp.Particle2Node(ts);

            sp.SolveSolid(ts);

            sp.Node2Particle(ts);

            sp.MoveParticle();

            if (ts.istep % params.iout == 0) {
                iview++;
                sp.OutputPointDataVTKHDF(iview, ts);
            }

            if (ts.istep % 100 == 0 && myrank == 0) { OutputMessage(iview, ts); }
        }

        if (params.rstflag == 2 || params.rstflag == 3) sp.RestartOutput();
    }

    FinalizeSimulation(start_time);

    return 0;
}
