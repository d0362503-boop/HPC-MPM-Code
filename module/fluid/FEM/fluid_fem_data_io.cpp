#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/shape_function.h"
#include "module/vtk_hdf5.h"

using namespace stabilizedfem;

void StabilizedFEM::InputBCData(std::ifstream &infile) {

    this->ubc.BCInput(infile);
    this->vbc.BCInput(infile);
    this->wbc.BCInput(infile);
    this->pbc.BCInput(infile);

    return;
}

void StabilizedFEM::RestartInput() {
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    this->InitializeMeshData();

    std::string filename = params.gridfile + std::to_string(myrank) + "_re.txt";

    std::ifstream reinfile;
    reinfile.open(filename);

    for (int n = 0; n < cmesh.nodec; n++) {
        reinfile >> this->naccel[n + cmesh.nuc] >> this->naccel[n + cmesh.nvc] >> this->naccel[n + cmesh.nwc]        //
            >> this->nvel[n + cmesh.nuc] >> this->nvel[n + cmesh.nvc] >> this->nvel[n + cmesh.nwc]                   //
            >> this->npres[n] >> this->nphi[n];
        reinfile.ignore(1000, '\n');
    }

    reinfile.close();

    return;
}

void StabilizedFEM::RestartOutput() {
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    std::string filename = params.gridfile + std::to_string(myrank) + "_re.txt";

    std::ofstream reoutfile;
    reoutfile.flags(std::ios::right | std::ios::scientific);
    reoutfile.open(filename);

    for (int n = 0; n < cmesh.nodec; n++) {
        reoutfile << std::setw(15) << this->naccel[n + cmesh.nuc]     //
                  << std::setw(15) << this->naccel[n + cmesh.nvc]     //
                  << std::setw(15) << this->naccel[n + cmesh.nwc]     //
                  << std::setw(15) << this->nvel[n + cmesh.nuc]       //
                  << std::setw(15) << this->nvel[n + cmesh.nvc]       //
                  << std::setw(15) << this->nvel[n + cmesh.nwc]       //
                  << std::setw(15) << this->npres[n]            //
                  << std::setw(15) << this->nphi[n] << "\n";
    }

    reoutfile.close();

    return;
}

void StabilizedFEM::OutputMeshDataVTKHDF(int iview, const TimeState &ts) {
    const Mesh &mesh = this->grid_.mesh;
    const RunParams &params = this->params_;

#ifdef HAVE_HDF5
    std::string filename = params.outfile + "-" + std::to_string(iview) + "-w.vtkhdf";

    std::vector<std::array<double, 3>> velocity(mesh.node);

    for (int n = 0; n < mesh.node; ++n) {
        velocity[n] = {this->nvel_vtk[n + mesh.nu], this->nvel_vtk[n + mesh.nv], this->nvel_vtk[n + mesh.nw]};
    }

    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_WORLD);

    // Compute global point offset so local connectivity is written as global node IDs.
    hsize_t local_npts = static_cast<hsize_t>(mesh.xyn.size());
    hsize_t point_global_offset = 0;
    hsize_t total_npts = local_npts;
    vtkhdf::VTKHDFWriter::ComputeGlobalInfo(local_npts, point_global_offset, total_npts, MPI_COMM_WORLD);

    auto info = vtkhdf::WriteHexMeshTopology(writer, mesh.xyn, mesh.nc, point_global_offset);
    writer.SetTime(ts.real_time);

    writer.CreatePointDataGroup();
    writer.WritePointVector("Velocity", info.total_npts, info.local_npts, info.point_global_offset, velocity);
    writer.WritePointScalar("Pressure", info.total_npts, info.local_npts, info.point_global_offset, this->npres_vtk);
    writer.WritePointScalar("Phi", info.total_npts, info.local_npts, info.point_global_offset, this->nphi_vtk);
#endif

    return;
}

void StabilizedFEM::Cp2NodeVTK() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(mesh.node, this->nphi_vtk);
    VectorAssign(mesh.node, this->npres_vtk);
    VectorAssign(mesh.node * 3, this->nvel_vtk);
    for (int m = 0; m < mesh.nelem; m++) {
        for (int n = 0; n < 8; n++) {
            int id = mesh.nc[m][n];
            std::array<double, 3> xyp = {mesh.xyn[id][0], mesh.xyn[id][1], mesh.xyn[id][2]};
            this->nphi_vtk[id] = 0.0e0;
            this->npres_vtk[id] = 0.0e0;
            this->nvel_vtk[id + mesh.nu] = 0.0e0;
            this->nvel_vtk[id + mesh.nv] = 0.0e0;
            this->nvel_vtk[id + mesh.nw] = 0.0e0;
            MakeSF(m, xyp, cmesh.idimc, cmesh.xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->nphi_vtk[id] += sfi * this->nphi[nid];
                this->npres_vtk[id] += sfi * this->npres[nid];
                this->nvel_vtk[id + mesh.nu] += sfi * this->nvel[nid + cmesh.nuc];
                this->nvel_vtk[id + mesh.nv] += sfi * this->nvel[nid + cmesh.nvc];
                this->nvel_vtk[id + mesh.nw] += sfi * this->nvel[nid + cmesh.nwc];
            }
        }
    }

    for (int n = 0; n < mesh.node; n++) {
        if (this->nphi_vtk[n] < 0.01e0) this->nphi_vtk[n] = 0.0e0;
        if (this->nphi_vtk[n] > 0.99e0) this->nphi_vtk[n] = 1.0e0;
    }

    return;
}
