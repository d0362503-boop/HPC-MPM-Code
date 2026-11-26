#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

using namespace stabilizedmpm;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        Grid grid;          // rank grid
        RunParams params;   // run parameters
        TimeState ts;       // time-step state

        StabilizedMPM wp(grid, params);

        wp.DataInput(grid, params);

        if (params.rstflag == 1 || params.rstflag == 3) wp.RestartInput();

        BuildMesh(grid.mesh);

        BuildControlPoint(grid.mesh, grid.cmesh);

        wp.NS_.BuildCrsMat(16);

        ComputeNodalVol(params, grid);

        ts.istep = params.ista - 1;
        int iview = ts.istep / params.iout;
        ts.real_time = params.dt * double(ts.istep);

        wp.OutputPointDataVTKHDF(iview, ts);

        for (ts.istep = params.ista; ts.istep <= params.iend; ts.istep++) {

            // -----------------------------------------------
            ts.real_time = params.dt * double(ts.istep);
            if (ts.istep <= params.nlstep) {
                ts.facl = params.dlstep * double(ts.istep);
            } else {
                ts.facl = 1.0e0;
            }
            // -----------------------------------------------

            wp.MeshPointLinklist();

            wp.Particle2Node(ts);

            wp.SolveNS(ts);

            wp.Node2Particle(ts);

            const bool is_dlb_step = wp.do_dlb && (ts.istep % params.iout == 0);

            (nprocs != 1 && is_dlb_step) ? wp.ApplyDLB(grid) : wp.MoveParticle();

            if (ts.istep % params.iout == 0) {
                iview++;
                wp.OutputPointDataVTKHDF(iview, ts);

                if (params.rstflag == 2 || params.rstflag == 3) wp.RestartOutput();
            }

            if (ts.istep % 100 == 0 && myrank == 0) { OutputMessage(iview, ts); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
