#pragma once

#ifdef HAVE_HDF5

#include <mpi.h>

#include <string>
#include <vector>

#include "material_point.h"
#include "vtk_hdf5.h"

// ----------------------------------------------------------------------------
// MPM-Code specific adapters: write MaterialPoint data as VTK HDF5.
// ----------------------------------------------------------------------------

namespace vtkhdf {

// Write all material points owned by this rank into a single collective
// VTK HDF5 UnstructuredGrid file.  Each particle becomes a VTK_VERTEX cell.
inline void OutputMaterialPointsToVtkHdf(const std::string& filename, const MaterialPoint& point, double time,
                                         MPI_Comm comm = MPI_COMM_WORLD) {
    int my_rank = 0;
    int n_ranks = 1;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &n_ranks);

    hsize_t local_npts = static_cast<hsize_t>(point.num);
    hsize_t global_offset = 0;
    hsize_t total_npts = 0;
    Writer::ComputeGlobalInfo(local_npts, global_offset, total_npts, comm);

    // Topology: each material point is a VTK_VERTEX cell.
    std::vector<long long> connectivity(local_npts);
    std::vector<long long> offsets(local_npts + 1);
    std::vector<unsigned char> types(local_npts, 1);  // VTK_VERTEX
    std::vector<std::array<double, 3>> points(local_npts);

    for (hsize_t i = 0; i < local_npts; ++i) {
        connectivity[i] = static_cast<long long>(global_offset + i);
        offsets[i] = static_cast<long long>(global_offset + i);
        points[i] = point.coord[i];
    }
    // Last offset value equals total number of points for vertex cells.
    offsets[local_npts] = static_cast<long long>(total_npts);

    Writer writer(filename, comm);
    writer.CreateUnstructuredGridGroup(total_npts, total_npts, total_npts);
    writer.SetTime(time);
    writer.WriteUnstructuredGridTopology(total_npts, local_npts, global_offset, points, connectivity, offsets, types);

    writer.CreatePointDataGroup();
    writer.WritePointScalarInt("matid", total_npts, local_npts, global_offset, point.matid);
    writer.WritePointVector("velocity", total_npts, local_npts, global_offset, point.vel);
    writer.WritePointVector("displacement", total_npts, local_npts, global_offset, point.displ);
    writer.WritePointVector("acceleration", total_npts, local_npts, global_offset, point.accel);
    writer.WritePointScalar("pressure", total_npts, local_npts, global_offset, point.pres);
    writer.WritePointScalar("volume", total_npts, local_npts, global_offset, point.vol);

    if (!point.stress.empty()) {
        std::vector<double> mean_stress(local_npts);
        std::vector<double> vm_stress(local_npts);
        for (hsize_t i = 0; i < local_npts; ++i) {
            const auto& s = point.stress[i];
            mean_stress[i] = (s[0] + s[1] + s[2]) / 3.0;
            vm_stress[i] = std::sqrt(0.5 * ((s[0] - s[1]) * (s[0] - s[1]) + (s[1] - s[2]) * (s[1] - s[2]) +
                                            (s[2] - s[0]) * (s[2] - s[0])) +
                                     3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
        }
        writer.WritePointScalar("mean_stress", total_npts, local_npts, global_offset, mean_stress);
        writer.WritePointScalar("von_mises_stress", total_npts, local_npts, global_offset, vm_stress);
    }

    writer.Close();
}

}  // namespace vtkhdf

#endif  // HAVE_HDF5
