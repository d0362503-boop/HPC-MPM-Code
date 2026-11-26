#include "data/generate/output_util.h"

#include <iomanip>
#include <string>

#include "module/mesh.h"

#ifdef HAVE_HDF5
#include <mpi.h>

#include "module/vtk_hdf5.h"
#endif

void WriteGlobalMeshHeader(std::ofstream &outfile, const Grid &grid) {
    const Mesh &mesh = grid.mesh;
    const ControlMesh &cmesh = grid.cmesh;

    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << cmesh.idimc[i];
    outfile << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << mesh.xymin[i];
    outfile << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << mesh.xymax[i];
    outfile << "\n";

    outfile << std::setw(10) << mesh.nelem;
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << mesh.xyelem[i];
    outfile << "\n";

    outfile << std::setw(10) << mesh.node;
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << mesh.xynode[i];
    outfile << "\n";

    outfile << std::setw(10) << cmesh.nodec;
    for (int i = 0; i < 3; ++i) outfile << std::setw(10) << cmesh.xynodec[i];
    outfile << "\n";
}

#ifdef HAVE_HDF5
void WriteVTKHDFMesh(const std::string &filename, const Mesh &mesh) {
    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_SELF);
    vtkhdf::WriteHexMeshTopology(writer, mesh.xyn, mesh.nc);
    writer.SetTime(0.0e0);
}

void WriteVTKHDFPoints(const std::string &filename, const MaterialPoint &point) {
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
