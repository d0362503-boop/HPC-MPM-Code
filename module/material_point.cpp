#include "material_point.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

#include "dataset.h"
#include "map_and_interpolate.h"
#include "mesh.h"
#include "mpi_data.h"
#include "shape_function.h"

void MaterialPoint::BuildGaussianPoint() {
    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    int npe = npxye[0] * npxye[1] * npxye[2];
    VectorAssign(nelem * npe, this->coord);

    int nex = xyelem[0], ney = xyelem[1], nez = xyelem[2];

    for (int m = 0; m < nelem; m++) {
        int ize = m / (nex * ney);
        int iye = (m - ize * (nex * ney)) / nex;
        int ixe = m - ize * (nex * ney) - iye * nex;

        std::array<double, 3> xye, xyp;
        xye[0] = xymin[0] + dxy[0] * (double(ixe) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(iye) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ize) + 0.5e0);
        for (int iz = 0; iz < npxye[2]; iz++) {
            xyp[2] = xye[2] + dec2p[iz][2];
            for (int iy = 0; iy < npxye[1]; iy++) {
                xyp[1] = xye[1] + dec2p[iy][1];
                for (int ix = 0; ix < npxye[0]; ix++) {
                    xyp[0] = xye[0] + dec2p[ix][0];
                    this->coord[this->num][0] = xyp[0];
                    this->coord[this->num][1] = xyp[1];
                    this->coord[this->num][2] = xyp[2];
                    this->num++;
                }
            }
        }
    }

    return;
}

void MaterialPoint::MeshPointLinklist() {
    std::vector<int> idepl(nelem, -1);
    VectorAssign(nelem, this->numep);
    VectorAssign(nelem, this->idepf, -1);
    VectorAssign(this->num, this->idp2p, -1);

    for (int ip = 0; ip < this->num; ip++) {
        double xx = this->coord[ip][0] - xymin[0];
        double yy = this->coord[ip][1] - xymin[1];
        double zz = this->coord[ip][2] - xymin[2];
        int iex = std::floor(xx / dxy[0]);
        int iey = std::floor(yy / dxy[1]);
        int iez = std::floor(zz / dxy[2]);
        int ne = iex + iey * xyelem[0] + iez * xyelem[0] * xyelem[1];
        this->numep[ne]++;
        if (this->numep[ne] == 1) {
            this->idepf[ne] = ip;
            idepl[ne] = ip;
        } else {
            int ipo = idepl[ne];
            this->idp2p[ipo] = ip;
            idepl[ne] = ip;
        }
    }

    return;
}

void MaterialPoint::CalPointUnitNormal() {
    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec * 3, this->nnormal);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = {this->coord[pid][0], this->coord[pid][1], this->coord[pid][2]};
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                this->nnormal[nid + nuc] += dsfi1 * this->mass[pid];
                this->nnormal[nid + nvc] += dsfi2 * this->mass[pid];
                this->nnormal[nid + nwc] += dsfi3 * this->mass[pid];
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nnormal, {nuc, nvc, nwc});

    for (int n = 0; n < nodec; n++) {
        double norm =
            std::sqrt(pow(this->nnormal[n + nuc], 2) + pow(this->nnormal[n + nvc], 2) + pow(this->nnormal[n + nwc], 2));
        if (norm > mtol) {
            this->nnormal[n + nuc] /= norm;
            this->nnormal[n + nvc] /= norm;
            this->nnormal[n + nwc] /= norm;
        } else {
            this->nnormal[n + nuc] = 0.0e0;
            this->nnormal[n + nvc] = 0.0e0;
            this->nnormal[n + nwc] = 0.0e0;
        }
    }

    // --- TO DO ---
    // --- Need unit normal BC setting ---

    return;
}

std::array<std::array<double, 3>, 3> MaterialPoint::ComputeDeltaDefGrad(const std::vector<int> &nc, int nenode,
                                                                        double af_coeff,
                                                                        const std::vector<std::array<double, 3>> &dsf) {
    std::array<std::array<double, 3>, 3> def_grad{};
    for (int ni = 0; ni < nenode; ni++) {
        int nid = nc[ni];
        double dsfi1 = dsf[ni][0];
        double dsfi2 = dsf[ni][1];
        double dsfi3 = dsf[ni][2];
        def_grad[0][0] += dsfi1 * af_coeff * this->ndispl[nid + nuc];
        def_grad[0][1] += dsfi2 * af_coeff * this->ndispl[nid + nuc];
        def_grad[0][2] += dsfi3 * af_coeff * this->ndispl[nid + nuc];
        def_grad[1][0] += dsfi1 * af_coeff * this->ndispl[nid + nvc];
        def_grad[1][1] += dsfi2 * af_coeff * this->ndispl[nid + nvc];
        def_grad[1][2] += dsfi3 * af_coeff * this->ndispl[nid + nvc];
        def_grad[2][0] += dsfi1 * af_coeff * this->ndispl[nid + nwc];
        def_grad[2][1] += dsfi2 * af_coeff * this->ndispl[nid + nwc];
        def_grad[2][2] += dsfi3 * af_coeff * this->ndispl[nid + nwc];
    }
    ++def_grad[0][0];
    ++def_grad[1][1];
    ++def_grad[2][2];

    return def_grad;
}

void MaterialPoint::ImplicitDsfCorr(const std::vector<int> &nc, int nenode, //
                                    std::vector<std::array<double, 3>> &dsf) {
    std::array<std::array<double, 3>, 3> def_grad{};
    def_grad = this->ComputeDeltaDefGrad(nc, nenode, this->alpha_f, dsf);

    std::array<std::array<double, 3>, 3> def_grad_inv;
    def_grad_inv = InvMat3(def_grad);

    std::array<double, 3> dsf_temp{};

    for (int ni = 0; ni < nenode; ni++) {
        for (int j = 0; j < 3; j++) {
            dsf_temp[j] = 0.0e0;
            for (int i = 0; i < 3; i++) { dsf_temp[j] += dsf[ni][i] * def_grad_inv[i][j]; }
        }
        for (int j = 0; j < 3; j++) { dsf[ni][j] = dsf_temp[j]; }
    }

    return;
}

std::vector<std::array<double, 3>> MaterialPoint::DeltaCorrectionParticleShifting() {
    std::vector<double> nei(nodec, 0.0e0);
    std::vector<std::array<double, 3>> delta_corr;
    VectorAssign(this->num, delta_corr);

    double eu_norm = 0.0e0;
    for (int i = 0; i < nodec; i++) {
        nei[i] = std::max(0.0e0, -nvol[i] + this->nvof[i]);
        eu_norm += (nei[i] * nei[i]) * dbc[i];
    }

    MPI_Allreduce(MPI_IN_PLACE, &eu_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;
    std::vector<std::array<double, 3>> geup(this->num);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                geup[pid][0] += dsfi1 * nei[nid];
                geup[pid][1] += dsfi2 * nei[nid];
                geup[pid][2] += dsfi3 * nei[nid];
            }
            pid = this->idp2p[pid];
        }
    }

    double geup_dot = 0.0e0;
    for (int n = 0; n < this->num; n++) {
        geup[n][0] *= 2.0e0 * this->vol[n];
        geup[n][1] *= 2.0e0 * this->vol[n];
        geup[n][2] *= 2.0e0 * this->vol[n];
        for (int i = 0; i < 3; i++) { geup_dot += geup[n][i] * geup[n][i]; }
    }
    MPI_Allreduce(MPI_IN_PLACE, &geup_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double b_0 = (geup_dot < mtol) ? 0.0e0 : (eu_norm / geup_dot);

    for (int n = 0; n < this->num; n++) {
        delta_corr[n][0] = -b_0 * geup[n][0];
        delta_corr[n][1] = -b_0 * geup[n][1];
        delta_corr[n][2] = -b_0 * geup[n][2];
    }

    return delta_corr;
}
