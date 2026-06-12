#include "mpi_data.h"
#include "dataset.h"
#include "material_point.h"
#include "mesh.h"
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <vector>

void MaterialPoint::DetermineParticleRank() {

    std::vector<double> dxyr(3);
    for (int i = 0; i < 3; i++) { dxyr[i] = xymax[i] - xymin[i] - dxy[i]; }

    std::vector<int> mxy(3), checkmp;
    checkmp.assign(this->num, -1);
    this->par_comm_.ofp_id.assign(this->num, 0);
    this->par_comm_.nsp.assign(isb, 0);
    this->par_comm_.nrp.assign(isb, 0);
    this->par_comm_.nnmp = 0;
    this->par_comm_.nrmp = 0;
    for (int ip = 0; ip < this->num; ip++) {
        bool skip_ip = false;
        for (int i = 0; i < 3; i++) {
            if (this->coord[ip][i] < xyminw[i] || this->coord[ip][i] > xymaxw[i]) {
                this->par_comm_.ofp_id[this->par_comm_.nrmp] = ip;
                this->par_comm_.nrmp++;
                skip_ip = true;
                break;
            } else if (this->coord[ip][i] < xymin[i] - dxyr[i] || this->coord[ip][i] > xymax[i] + dxyr[i]) {
                this->par_comm_.nrmp++;
                std::cout << "movepts error" << "\n";
                exit(1);
            } else if (this->coord[ip][i] < xymin[i]) {
                mxy[i] = -1;
            } else if (this->coord[ip][i] > xymax[i]) {
                mxy[i] = 1;
            } else {
                mxy[i] = 0;
            }
        }

        if (skip_ip) continue;

        int idsb = myrank + mxy[0] + mxy[1] * nxyr[0] + mxy[2] * nxyr[0] * nxyr[1];

        if (idsb == myrank) {
            this->par_comm_.nnmp++;
            checkmp[ip] = -2;
        } else {
            for (int i = 0; i < isb; i++) {
                if (naid[i] == idsb) {
                    this->par_comm_.nsp[i]++;
                    checkmp[ip] = i;
                    break;
                }
            }
        }
    }

    std::vector<int> idmploc(isb);
    this->par_comm_.nmps = 0;
    for (int i = 0; i < isb; i++) {
        idmploc[i] = this->par_comm_.nmps;
        this->par_comm_.nmps += this->par_comm_.nsp[i];
    }

    VectorAssign(this->par_comm_.nmps, this->par_comm_.idmp);
    VectorAssign(this->par_comm_.nnmp, this->par_comm_.idnsp);

    int nnmp2 = 0;
    for (int ip = 0; ip < this->num; ip++) {
        int i = checkmp[ip];
        if (i == -2) { // Stay
            this->par_comm_.idnsp[nnmp2++] = ip;
        } else if (i == -1) { // Remove
            continue;
        } else { // Move
            this->par_comm_.idmp[idmploc[i]++] = ip;
        }
    }

    std::vector<MPI_Request> irqs(isb);
    std::vector<MPI_Request> irqr(isb);
    MPI_Status status;
    for (int i = 0; i < isb; i++) {
        int ncomid = naid[i];
        MPI_Isend(&this->par_comm_.nsp[i], 1, MPI_INT, ncomid, 1, MPI_COMM_WORLD, &irqs[i]);
        MPI_Irecv(&this->par_comm_.nrp[i], 1, MPI_INT, ncomid, 1, MPI_COMM_WORLD, &irqr[i]);
    }
    for (int i = 0; i < isb; i++) {
        MPI_Wait(&irqs[i], &status);
        MPI_Wait(&irqr[i], &status);
    }

    this->par_comm_.nrps = 0;
    for (int i = 0; i < isb; i++) { this->par_comm_.nrps += this->par_comm_.nrp[i]; }

    return;
}
