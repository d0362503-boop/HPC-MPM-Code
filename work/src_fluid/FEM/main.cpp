#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

using namespace stabilizedfem;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        StabilizedFEM fluid;

        fluid.DataInput();

        if (rstflag == 1 || rstflag == 3) fluid.RestartInput();

        BuildMesh(g_grid.mesh);

        BuildControlPoint(g_grid.mesh, g_grid.cmesh);

        fluid.BuildSubcellPoints();

        fluid.MeshPointLinklist();

        fluid.Particle2NodePhi();

        fluid.NS_.BuildCrsMat(16);

        ComputeNodalVol(g_params, g_grid);

        istep = ista - 1;
        int iview = istep / iout;
        real_time = dt * double(istep);

        fluid.Cp2NodeVTK();
        fluid.OutputMeshDataVTKHDF(iview, g_time);

        for (istep = ista; istep <= iend; istep++) {

            // -----------------------------------------------
            real_time = dt * double(istep);
            if (istep <= nlstep) {
                facl = dlstep * double(istep);
            } else {
                facl = 1.0e0;
            }
            // -----------------------------------------------

            fluid.SetPFDomain();

            fluid.SolveNS(g_time);

            fluid.CommitNodalKinematics();

            if (istep % iout == 0) {
                iview++;
                fluid.Cp2NodeVTK();
                fluid.OutputMeshDataVTKHDF(iview, g_time);

                if (rstflag == 2 || rstflag == 3) fluid.RestartOutput();
            }

            if (istep % 100 == 0 && myrank == 0) { OutputMessage(iview, g_time); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
