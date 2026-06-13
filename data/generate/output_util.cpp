#include "output_util.h"

#include <iomanip>
#include <string>

#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/mesh.h"

#ifdef HAVE_HDF5
#include <mpi.h>

#include "../../module/vtk_hdf5.h"
#endif

void WriteGlobalMeshHeader(std::ofstream &outfile) {
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << idimc[i];
    outfile << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << xymin[i];
    outfile << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << xymax[i];
    outfile << "\n";

    outfile << std::setw(10) << nelem;
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << xyelem[i];
    outfile << "\n";

    outfile << std::setw(10) << node;
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << xynode[i];
    outfile << "\n";

    outfile << std::setw(10) << nodec;
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << xynodec[i];
    outfile << "\n";
}

#ifdef HAVE_HDF5
void WriteVtkHdf5Mesh(const std::string &filename) {
    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_SELF);
    vtkhdf::WriteHexMeshTopology(writer, xyn, nc);
    writer.SetTime(0.0e0);
}

void WriteVtkHdf5Points(const std::string &filename, const MaterialPoint &point) {
    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_SELF);
    auto info = vtkhdf::WriteParticleTopology(writer, point.coord);
    writer.SetTime(0.0e0);

    writer.CreatePointDataGroup();
    if (!point.id.empty()) {
        writer.WritePointScalar("ID", info.total_npts, info.local_npts, info.global_offset, point.id);
    }
    if (!point.matid.empty()) {
        writer.WritePointScalar("MatID", info.total_npts, info.local_npts, info.global_offset, point.matid);
    }
    if (!point.surf_point.empty()) {
        writer.WritePointScalar("SurfFlag", info.total_npts, info.local_npts, info.global_offset, point.surf_point);
    }
}
#endif
