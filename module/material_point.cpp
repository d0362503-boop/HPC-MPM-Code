#include "module/material_point.h"

#include <mpi.h>

#include <cmath>
#include <vector>

#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

void MaterialPoint::BuildSubcellPoints() {
    const Mesh &mesh = this->grid_.mesh;
    const RunParams &params = this->params_;

    std::array<std::array<double, 3>, 6> dec2p;
    SubcellOffsets(params.npxye, mesh.dxy, dec2p);

    int npe = params.npxye[0] * params.npxye[1] * params.npxye[2];
    VectorAssign(mesh.nelem * npe, this->coord);

    for (int m = 0; m < mesh.nelem; m++) {
        const std::array<double, 3> xye = ElementCenter(mesh, m);
        std::array<double, 3> xyp;
        for (int iz = 0; iz < params.npxye[2]; iz++) {
            xyp[2] = xye[2] + dec2p[iz][2];
            for (int iy = 0; iy < params.npxye[1]; iy++) {
                xyp[1] = xye[1] + dec2p[iy][1];
                for (int ix = 0; ix < params.npxye[0]; ix++) {
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
    const Mesh &mesh = this->grid_.mesh;

    std::vector<int> idepl(mesh.nelem, -1);
    VectorAssign(mesh.nelem, this->numep);
    VectorAssign(mesh.nelem, this->idepf, -1);
    VectorAssign(this->num, this->idp2p, -1);

    for (int ip = 0; ip < this->num; ip++) {
        double xx = this->coord[ip][0] - mesh.xymin[0];
        double yy = this->coord[ip][1] - mesh.xymin[1];
        double zz = this->coord[ip][2] - mesh.xymin[2];
        int iex = static_cast<int>(std::floor(xx / mesh.dxy[0]));
        int iey = static_cast<int>(std::floor(yy / mesh.dxy[1]));
        int iez = static_cast<int>(std::floor(zz / mesh.dxy[2]));
        int ne = iex + iey * mesh.xyelem[0] + iez * mesh.xyelem[0] * mesh.xyelem[1];
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

void MaterialPoint::CalNodalUnitNormal() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(cmesh.nodec * 3, this->nnormal);
    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(this->grid_, m, xyp, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                this->nnormal[nid + cmesh.nuc] += dsfi1 * this->mass[pid];
                this->nnormal[nid + cmesh.nvc] += dsfi2 * this->mass[pid];
                this->nnormal[nid + cmesh.nwc] += dsfi3 * this->mass[pid];
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->grid_.part, this->nnormal, {cmesh.nuc, cmesh.nvc, cmesh.nwc});

    for (int n = 0; n < cmesh.nodec; n++) {
        double norm = std::sqrt(std::pow(this->nnormal[n + cmesh.nuc], 2)   //
                                + std::pow(this->nnormal[n + cmesh.nvc], 2) //
                                + std::pow(this->nnormal[n + cmesh.nwc], 2));
        if (norm > params.mtol) {
            this->nnormal[n + cmesh.nuc] /= norm;
            this->nnormal[n + cmesh.nvc] /= norm;
            this->nnormal[n + cmesh.nwc] /= norm;
        } else {
            this->nnormal[n + cmesh.nuc] = 0.0e0;
            this->nnormal[n + cmesh.nvc] = 0.0e0;
            this->nnormal[n + cmesh.nwc] = 0.0e0;
        }
    }

    return;
}

void MaterialPoint::CalPointUnitNormal(std::vector<std::array<double, 3>> &particle_normal, std::vector<double> &surface_weight) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->num, particle_normal);
    VectorAssign(this->num, surface_weight);
    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            MakeSF(this->grid_, m, this->coord[pid], ncm, nenode, sf, dsf);

            for (int ni = 0; ni < nenode; ni++) {
                const int nid = ncm[ni];
                const double sfi = sf[ni];
                const double phi = this->nphi[nid];
                surface_weight[pid] += sfi * phi;
                particle_normal[pid][0] += sfi * this->nnormal[nid + cmesh.nuc];
                particle_normal[pid][1] += sfi * this->nnormal[nid + cmesh.nvc];
                particle_normal[pid][2] += sfi * this->nnormal[nid + cmesh.nwc];
            }

            const std::array<double, 3> normal = particle_normal[pid];
            double normal_norm = NormVec3(normal);
            for (int i = 0; i < 3; i++) { particle_normal[pid][i] /= normal_norm; }

            pid = this->idp2p[pid];
        }
    }

    return;
}

std::array<std::array<double, 3>, 3>
MaterialPoint::ComputeDeltaDefGrad(const std::vector<int> &nc, int nenode, //
                                   double af_coeff, const std::vector<std::array<double, 3>> &dsf) const noexcept {
    const ControlMesh &cmesh = this->grid_.cmesh;

    std::array<std::array<double, 3>, 3> ddg{};
    for (int ni = 0; ni < nenode; ni++) {
        int nid = nc[ni];
        double dsfi1 = dsf[ni][0];
        double dsfi2 = dsf[ni][1];
        double dsfi3 = dsf[ni][2];
        ddg[0][0] += dsfi1 * af_coeff * this->ndispl[nid + cmesh.nuc];
        ddg[0][1] += dsfi2 * af_coeff * this->ndispl[nid + cmesh.nuc];
        ddg[0][2] += dsfi3 * af_coeff * this->ndispl[nid + cmesh.nuc];
        ddg[1][0] += dsfi1 * af_coeff * this->ndispl[nid + cmesh.nvc];
        ddg[1][1] += dsfi2 * af_coeff * this->ndispl[nid + cmesh.nvc];
        ddg[1][2] += dsfi3 * af_coeff * this->ndispl[nid + cmesh.nvc];
        ddg[2][0] += dsfi1 * af_coeff * this->ndispl[nid + cmesh.nwc];
        ddg[2][1] += dsfi2 * af_coeff * this->ndispl[nid + cmesh.nwc];
        ddg[2][2] += dsfi3 * af_coeff * this->ndispl[nid + cmesh.nwc];
    }

    for (int i = 0; i < 3; i++) { ++ddg[i][i]; }

    return ddg;
}

