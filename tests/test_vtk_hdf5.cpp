// Smoke test for the VTK HDF5 writer class.
// Writes a tiny UnstructuredGrid with 4 vertex cells per rank directly using
// the Writer class, without instantiating MaterialPoint.

#ifdef HAVE_HDF5

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <vector>

#include "vtk_hdf5.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int nprocs = 1;
    int myrank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    try {
        hsize_t local_n = 4;
        hsize_t global_offset = 0;
        hsize_t total_n = 0;
        vtkhdf::Writer::ComputeGlobalInfo(local_n, global_offset, total_n, MPI_COMM_WORLD);

        std::vector<long long> connectivity(local_n);
        std::vector<long long> offsets(local_n + 1);
        std::vector<unsigned char> types(local_n, 1);  // VTK_VERTEX
        std::vector<std::array<double, 3>> points(local_n);
        std::vector<double> scalar_field(local_n);
        std::vector<std::array<double, 3>> vector_field(local_n);

        for (hsize_t i = 0; i < local_n; ++i) {
            hsize_t idx = global_offset + i;
            connectivity[i] = static_cast<long long>(idx);
            offsets[i] = static_cast<long long>(idx);
            points[i] = {static_cast<double>(idx), static_cast<double>(idx) * 0.5,
                         static_cast<double>(idx) * 0.25};
            scalar_field[i] = static_cast<double>(myrank * 10 + i);
            vector_field[i] = {1.0, 0.0, static_cast<double>(idx)};
        }
        offsets[local_n] = static_cast<long long>(total_n);

        std::string filename = "test_grid.vtkhdf";
        vtkhdf::Writer writer(filename, MPI_COMM_WORLD);
        writer.CreateUnstructuredGridGroup(total_n, total_n, total_n);
        writer.SetTime(0.0);
        writer.WriteUnstructuredGridTopology(total_n, local_n, global_offset, points, connectivity, offsets, types);

        writer.CreatePointDataGroup();
        writer.WritePointScalar("scalar", total_n, local_n, global_offset, scalar_field);
        writer.WritePointVector("vector", total_n, local_n, global_offset, vector_field);
        writer.Close();

        if (myrank == 0) {
            std::cout << "VTK HDF5 smoke test passed: total points = " << total_n << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "[rank " << myrank << "] VTK HDF5 test failed: " << e.what() << std::endl;
        MPI_Finalize();
        return 1;
    }

    MPI_Finalize();
    return 0;
}

#else

int main() {
    std::cout << "HDF5 support disabled; test skipped." << std::endl;
    return 0;
}

#endif
