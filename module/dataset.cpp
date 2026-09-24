#include "module/dataset.h"

#include <iomanip>
#include <iostream>
#include <petsc.h>

#include "module/mpi_data.h"

double InitializeSimulation(int &argc, char **&argv) {

    PetscInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);

    return MPI_Wtime();
}

void FinalizeSimulation(double start_time) {

    const double calc_time = MPI_Wtime() - start_time;
    if (myrank == 0) {
        std::cout << "Calculation time: " << std::setw(15) << calc_time << "\n";
        std::cout << " ---- Job Finished ----" << "\n";
    }
    PetscFinalize();

    return;
}
