#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

using namespace stabilizedfem;

void StabilizedFEM::Particle2NodePhi() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    VectorAssign(this->num, this->phi);
    for (int ip = 0; ip < this->num; ip++) { this->phi[ip] = (this->coord[ip][2] < fs_height) ? 1.0e0 : 0.0e0; }

    double volp = (mesh.dxy[0] * mesh.dxy[1] * mesh.dxy[2]) / (params.npxye[0] * params.npxye[1] * params.npxye[2]);

    VectorAssign(cmesh.nodec, this->nvof);
    VectorAssign(cmesh.nodec, this->nphi);
    for (int m = 0; m < mesh.nelem; m++) {

        int nenode;
        std::vector<int> ncm;
        std::vector<double> sf;
        std::vector<std::array<double, 3>> dsf;

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = {this->coord[pid][0], this->coord[pid][1], this->coord[pid][2]};
            MakeSF(this->grid_, m, xyp, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->nvof[nid] += sfi * volp;
                this->nphi[nid] += sfi * this->phi[pid] * volp;
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nvof, 0);
    NodeVarComm(this->nphi, 0);

    for (int n = 0; n < cmesh.nodec; n++) {
        this->nphi[n] /= this->nvof[n];
        if (this->nphi[n] < 0.0e0) {
            this->nphi[n] = 0.0e0;
        } else if (this->nphi[n] > 1.0e0) {
            this->nphi[n] = 1.0e0;
        }
    }

    return;
}

double StabilizedFEM::CalLiquidVol() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    int nenode = (cmesh.idimc[0] + 1) * (cmesh.idimc[1] + 1) * (cmesh.idimc[2] + 1);

    double vol_liquid = 0.0e0;
    for (int m = 0; m < mesh.nelem; m++) {
        double mesh_vol = mesh.dxy[0] * mesh.dxy[1] * mesh.dxy[2];

        double mesh_phi = 0.0e0;
        for (int n = 0; n < nenode; n++) {
            int nid = cmesh.ncc[m][n];
            mesh_phi += this->nphi[nid] / double(nenode);
        }
        vol_liquid += mesh_phi * mesh_vol;
    }

    MPI_Allreduce(MPI_IN_PLACE, &vol_liquid, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return vol_liquid;
}

void StabilizedFEM::SetPFDomain() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    for (int n = 0; n < cmesh.nodec; n++) {
        if (this->nphi[n] < 0.0e0) {
            this->nphi[n] = 0.0e0;
        } else if (this->nphi[n] > 1.0e0) {
            this->nphi[n] = 1.0e0;
        }
    }

    int nenode = (cmesh.idimc[0] + 1) * (cmesh.idimc[1] + 1) * (cmesh.idimc[2] + 1);

    VectorAssign(mesh.nelem, this->rhoe);
    VectorAssign(mesh.nelem, this->rmue);
    for (int m = 0; m < mesh.nelem; m++) {
        for (int n = 0; n < nenode; n++) {
            int nid = cmesh.ncc[m][n];
            this->rhoe[m] += (this->nphi[nid] * this->rhol + (1.0e0 - this->nphi[nid]) * this->rhog) / double(nenode);
            this->rmue[m] += (this->nphi[nid] * this->rmul + (1.0e0 - this->nphi[nid]) * this->rmug) / double(nenode);
        }
    }

    return;
}
