#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/explicit/explicit_mpm_solid.h"

using namespace explicitmpm;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        ExplicitSolidMPM sp;

        sp.DataInput();

        if (rstflag == 1 || rstflag == 3) sp.RestartInput();

        BuildMesh(g_grid.mesh);

        BuildControlPoint(g_grid.mesh, g_grid.cmesh);

        ComputeNodalVol(g_params, g_grid);

        istep = ista - 1;
        int iview = istep / iout;
        real_time = dt * double(istep);

        sp.OutputPointDataVTKHDF(iview, g_time);

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

            sp.Particle2Node(g_time);

            sp.SolveSolid(g_time);

            sp.Node2Particle(g_time);

            sp.MoveParticle();

            if (istep % iout == 0) {
                iview++;
                sp.OutputPointDataVTKHDF(iview, g_time);
            }

            if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, g_time); }
        }

        if (rstflag == 2 || rstflag == 3) sp.RestartOutput();
    }

    FinalizeSimulation(start_time);

    return 0;
}
