#include "module/material_point.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

#include "module/dataset.h"
#include "module/map_and_interpolate.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

void MaterialPoint::BuildGaussianPoint() {
    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    int npe = npxye[0] * npxye[1] * npxye[2];
    VectorAssign(nelem * npe, this->coord);

    for (int m = 0; m < nelem; m++) {
        const std::array<int, 3> ijk = IndexToIJK(m, xyelem);

        std::array<double, 3> xye, xyp;
        xye[0] = xymin[0] + dxy[0] * (double(ijk[0]) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(ijk[1]) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ijk[2]) + 0.5e0);
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
        int iex = static_cast<int>(std::floor(xx / dxy[0]));
        int iey = static_cast<int>(std::floor(yy / dxy[1]));
        int iez = static_cast<int>(std::floor(zz / dxy[2]));
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
            std::array<double, 3> xyp = this->coord[pid];
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
        double norm = std::sqrt(std::pow(this->nnormal[n + nuc], 2)   //
                                + std::pow(this->nnormal[n + nvc], 2) //
                                + std::pow(this->nnormal[n + nwc], 2));
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

std::array<std::array<double, 3>, 3>
MaterialPoint::ComputeDeltaDefGrad(const std::vector<int> &nc, int nenode, double af_coeff,
                                   const std::vector<std::array<double, 3>> &dsf) const noexcept {

    std::array<std::array<double, 3>, 3> ddg{};
    for (int ni = 0; ni < nenode; ni++) {
        int nid = nc[ni];
        double dsfi1 = dsf[ni][0];
        double dsfi2 = dsf[ni][1];
        double dsfi3 = dsf[ni][2];
        ddg[0][0] += dsfi1 * af_coeff * this->ndispl[nid + nuc];
        ddg[0][1] += dsfi2 * af_coeff * this->ndispl[nid + nuc];
        ddg[0][2] += dsfi3 * af_coeff * this->ndispl[nid + nuc];
        ddg[1][0] += dsfi1 * af_coeff * this->ndispl[nid + nvc];
        ddg[1][1] += dsfi2 * af_coeff * this->ndispl[nid + nvc];
        ddg[1][2] += dsfi3 * af_coeff * this->ndispl[nid + nvc];
        ddg[2][0] += dsfi1 * af_coeff * this->ndispl[nid + nwc];
        ddg[2][1] += dsfi2 * af_coeff * this->ndispl[nid + nwc];
        ddg[2][2] += dsfi3 * af_coeff * this->ndispl[nid + nwc];
    }
    ++ddg[0][0];
    ++ddg[1][1];
    ++ddg[2][2];

    return ddg;
}

void MaterialPoint::ImplicitDsfCorr(const std::vector<int> &nc, int nenode, //
                                    std::vector<std::array<double, 3>> &dsf) const noexcept {
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
        dsf[ni] = dsf_temp;
    }

    return;
}
