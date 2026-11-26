#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

using namespace stabilizedfem;

int main(int argc, char *argv[]) {

    const double start_time = InitializeSimulation(argc, argv);

    {
        // Transitional: CrsMat still reads the global grid, run parameters and step index, so the
        // driver works on the global instances until CrsMat takes them explicitly.
        Grid &grid = g_grid;
        RunParams &params = g_params;
        TimeState &ts = g_time;

        StabilizedFEM fluid(grid, params);

        fluid.DataInput(grid, params);

        if (params.rstflag == 1 || params.rstflag == 3) fluid.RestartInput();

        BuildMesh(grid.mesh);

        BuildControlPoint(grid.mesh, grid.cmesh);

        fluid.BuildSubcellPoints();

        fluid.MeshPointLinklist();

        fluid.Particle2NodePhi();

        fluid.NS_.BuildCrsMat(16);

        ComputeNodalVol(params, grid);

        ts.istep = params.ista - 1;
        int iview = ts.istep / params.iout;
        ts.real_time = params.dt * double(ts.istep);

        fluid.Cp2NodeVTK();
        fluid.OutputMeshDataVTKHDF(iview, ts);

        for (ts.istep = params.ista; ts.istep <= params.iend; ts.istep++) {

            // -----------------------------------------------
            ts.real_time = params.dt * double(ts.istep);
            if (ts.istep <= params.nlstep) {
                ts.facl = params.dlstep * double(ts.istep);
            } else {
                ts.facl = 1.0e0;
            }
            // -----------------------------------------------

            fluid.SetPFDomain();

            fluid.SolveNS(ts);

            fluid.CommitNodalSolution();

            if (ts.istep % params.iout == 0) {
                iview++;
                fluid.Cp2NodeVTK();
                fluid.OutputMeshDataVTKHDF(iview, ts);

                if (params.rstflag == 2 || params.rstflag == 3) fluid.RestartOutput();
            }

            if (ts.istep % 100 == 0 && myrank == 0) { OutputMessage(iview, ts); }
        }
    }

    FinalizeSimulation(start_time);

    return 0;
}
