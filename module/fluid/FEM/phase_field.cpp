#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <optional>
#include <string>
#include <mpi.h>
#include <cmath>
#include "../../dataset.h"
#include "../../shape_function.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "stabilized_fem.h"

using namespace stabilizedfem;

void StabilizedFEM::Particle2NodePhi() {

    VectorAssign(this->num, this->phi);
    for (int ip = 0; ip < this->num; ip++) {
        this->phi[ip] = (this->coord[ip][2] < fs_height) ? 1.0e0 : 0.0e0;
    }

    double volp = (dxy[0] * dxy[1] * dxy[2]) / (npxye[0] * npxye[1] * npxye[2]);

    VectorAssign(nodec, this->nvof);
    VectorAssign(nodec, this->nphi);
    for (int m = 0; m < nelem; m++) {

        int nenode;
        std::vector<int> ncm; 
        std::vector<double> sf;
        std::vector<std::array<double, 3>> dsf;

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = {this->coord[pid][0], this->coord[pid][1], this->coord[pid][2]};
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
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

    for (int n = 0; n < nodec; n++) {
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

    int nenode = (idimc[0] + 1) * (idimc[1] + 1) * (idimc[2] + 1);

    double vol_liquid = 0.0e0;
    for (int m = 0; m < nelem; m++) {
        double mesh_vol = dxy[0] * dxy[1] * dxy[2];

        double mesh_phi = 0.0e0;
        for (int n = 0; n < nenode; n++) {
            int nid = ncc[m][n];
            mesh_phi += this->nphi[nid] / double(nenode);
        }
        vol_liquid += mesh_phi * mesh_vol;
    }

    MPI_Allreduce(MPI_IN_PLACE, &vol_liquid, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return vol_liquid;
}

void StabilizedFEM::SetPFDomain () {

    for (int n = 0; n < nodec; n++) {
        if (this->nphi[n] < 0.0e0) {
            this->nphi[n] = 0.0e0;
        } else if (this->nphi[n] > 1.0e0) {
            this->nphi[n] = 1.0e0;
        }
    }

    int nenode = (idimc[0] + 1) * (idimc[1] + 1) * (idimc[2] + 1);

    VectorAssign(nelem, this->rhoe);
    VectorAssign(nelem, this->rmue);
    for (int m = 0; m < nelem; m++) {
        for (int n = 0; n < nenode; n++) {
            int nid = ncc[m][n];
            this->rhoe[m] += (this->nphi[nid] * this->rhol 
                          +  (1.0e0 - this->nphi[nid]) * this->rhog) / double(nenode); 
            this->rmue[m] += (this->nphi[nid] * this->rmul 
                          +  (1.0e0 - this->nphi[nid]) * this->rmug) / double(nenode); 
        }
    }

    return;
}
