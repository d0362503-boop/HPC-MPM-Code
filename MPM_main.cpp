#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/mpi_data.h"
#include "module/solid/explicit/explicit_mpm_solid.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "work/src_fsi/MPM_FEM/block_fsi.h"
#include "work/src_fsi/MPM_MPM/block_fsi.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <petsc.h>

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
    double start_time = MPI_Wtime();

    // stabilizedmpm::StabilizedMixedMPM();

    // explicitmpm::SolidExplicitULMPM();

    // implicitmpm::SolidImplicitULMPM();

    // mpmfemblockfsi::ImmersedMPMFEMBlockFSI();

    mpmmpmblockfsi::MPMBlockFSI();

    double end_time = MPI_Wtime();
    double calc_time = end_time - start_time;
    if (myrank == 0) {
        std::cout << "Calculation time: " << std::setw(15) << calc_time << "\n";
        std::cout << " ---- Job Finished ----" << "\n";
    }
    PetscFinalize();

    return 0;
}
