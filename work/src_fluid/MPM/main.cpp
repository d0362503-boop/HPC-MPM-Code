#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

using namespace stabilizedmpm;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        StabilizedMPM wp;

        wp.DataInput();

        if (rstflag == 1 || rstflag == 3) wp.RestartInput();

        BuildMesh(g_grid.mesh);

        BuildControlPoint(g_grid.mesh, g_grid.cmesh);

        wp.NS_.BuildCrsMat(16);

        ComputeNodalVol(g_params, g_grid);

        istep = ista - 1;
        int iview = istep / iout;
        real_time = dt * double(istep);

        wp.OutputPointDataVTKHDF(iview, g_time);

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

            wp.Particle2Node(g_time);

            wp.SolveNS(g_time);

            wp.Node2Particle(g_time);

            const bool is_dlb_step = wp.do_dlb && (istep % iout == 0);

            (nprocs != 1 && is_dlb_step) ? wp.ApplyDLB() : wp.MoveParticle();

            if (istep % iout == 0) {
                iview++;
                wp.OutputPointDataVTKHDF(iview, g_time);

                if (rstflag == 2 || rstflag == 3) wp.RestartOutput();
            }

            if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, g_time); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
