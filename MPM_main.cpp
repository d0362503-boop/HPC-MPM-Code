#include "module/dataset.h"
#include "module/mpi_data.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <petsc.h>

// -----------------------
void StabilizedMixedMPM();
void Solid_Explicit_ULMPM();
void Solid_implicit_ULMPM();
void MPMBlockFSI();
// -----------------------

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
    double start_time = MPI_Wtime();

    // StabilizedMixedMPM();

    // Solid_Explicit_ULMPM();

    // Solid_implicit_ULMPM();  // ---- currently only for linear elasticity ----

    MPMBlockFSI(); // --- MPM FSI strong coupling with block iterative ---

    double end_time = MPI_Wtime();
    double calc_time = end_time - start_time;
    if (myrank == 0) {
        std::cout << "Calculation time: " << std::setw(15) << calc_time << "\n";
        std::cout << " ---- Job Finished ----" << "\n";
    }
    PetscFinalize();

    return 0;
}
