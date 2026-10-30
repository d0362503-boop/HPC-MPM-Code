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

void MaterialPoint::CommitNodalKinematics(const std::vector<double> &nvel_k, const std::vector<double> &naccel_k) {

    std::copy(nvel_k.begin(), nvel_k.end(), this->nvel.begin());
    std::copy(naccel_k.begin(), naccel_k.end(), this->naccel.begin());

    return;
}

std::array<bool, 6> MaterialPoint::DetermineClosedBoundaryFaces() const {

    // Count stationary normal BCs on each global face.
    std::array<double, 6> zero_normal_nodes{};
    std::array<int, 6> inflow_faces{};
    for (int i = 0; i < 3; i++) {
        const BoundaryCondition &normal_bc = i == 0 ? this->ubc : (i == 1 ? this->vbc : this->wbc);
        const BoundaryCondition &inflow_bc = i == 0 ? this->uinfbc : (i == 1 ? this->vinfbc : this->winfbc);
        const int lower_face = 2 * i;
        const int upper_face = lower_face + 1;
        for (int k = 0; k < normal_bc.ibc; k++) {
            if (normal_bc.fbc[k] != 0.0) continue;
            const int nid = normal_bc.nbc[k];
            const int index = IndexToIJK(nid, xynodec)[i] + aelemmin[i];
            if (index == 0) zero_normal_nodes[lower_face] += dbc[nid];
            if (index == xynodecw[i] - 1) zero_normal_nodes[upper_face] += dbc[nid];
        }
        for (int cell : inflow_bc.nbc) {
            const int index = IndexToIJK(cell, xyelem)[i] + aelemmin[i];
            if (index == 0) inflow_faces[lower_face] = 1;
            if (index == xyelemw[i] - 1) inflow_faces[upper_face] = 1;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, zero_normal_nodes.data(), 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, inflow_faces.data(), 6, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Close only fully constrained faces without normal inflow.
    std::array<bool, 6> closed_face{};
    for (int i = 0; i < 3; i++) {
        const int face_nodes = xynodecw[(i + 1) % 3] * xynodecw[(i + 2) % 3];
        for (int side = 0; side < 2; side++) {
            const int face = 2 * i + side;
            closed_face[face] = std::llround(zero_normal_nodes[face]) == face_nodes && !inflow_faces[face];
        }
    }

    return closed_face;
}

double MaterialPoint::ConstrainClosedBoundaryCoordinate(double coordinate, int direction,
                                                        const std::array<bool, 6> &closed_face) const {

    if (coordinate >= xyminw[direction] && coordinate <= xymaxw[direction]) return coordinate;
    const bool below_wall = coordinate < xyminw[direction];
    const int face = 2 * direction + (below_wall ? 0 : 1);
    if (!closed_face[face]) return coordinate;
    const double sign = below_wall ? 1.0 : -1.0;
    const double wall = below_wall ? xyminw[direction] : xymaxw[direction];

    return wall + sign * 0.01e0 * dxy[direction];
}

void MaterialPoint::CommitImplicitParticleKinematics(const std::vector<std::array<double, 3>> &accel_old,
                                                     const std::vector<std::array<double, 3>> &disp,
                                                     const std::vector<std::array<double, 3>> &disp_corr) {

    // const auto closed_face = this->DetermineClosedBoundaryFaces();

    const bool has_shift = !disp_corr.empty();
    for (int n = 0; n < this->num; n++) {
        for (int i = 0; i < 3; i++) {
            if (this->solswitch == MapScheme::FLIP) {
                this->vel[n][i] += dt * ((1.0e0 - this->integrator_.gamma_nb) * accel_old[n][i] //
                                         + this->integrator_.gamma_nb * this->accel[n][i]);
            }
            this->coord[n][i] += disp[n][i];
            if (has_shift) this->coord[n][i] += disp_corr[n][i];
            // this->coord[n][i] = this->ConstrainClosedBoundaryCoordinate(this->coord[n][i], i, closed_face);
        }
    }

    return;
}

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

void MaterialPoint::CalNodalUnitNormal() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec * 3, this->nnormal);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
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

    return;
}

void MaterialPoint::CalPointUnitNormal(std::vector<std::array<double, 3>> &particle_normal, std::vector<double> &surface_weight) {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->num, particle_normal);
    VectorAssign(this->num, surface_weight);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            MakeSF(m, this->coord[pid], idimc, xynodec, ncm, nenode, sf, dsf);

            for (int ni = 0; ni < nenode; ni++) {
                const int nid = ncm[ni];
                const double sfi = sf[ni];
                const double phi = this->nphi[nid];
                surface_weight[pid] += sfi * phi;
                particle_normal[pid][0] += sfi * this->nnormal[nid + nuc];
                particle_normal[pid][1] += sfi * this->nnormal[nid + nvc];
                particle_normal[pid][2] += sfi * this->nnormal[nid + nwc];
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

    for (int i = 0; i < 3; i++) { ++ddg[i][i]; }

    return ddg;
}

void MaterialPoint::ImplicitDsfCorr(const std::vector<int> &nc, int nenode, //
                                    std::vector<std::array<double, 3>> &dsf) const noexcept {

    std::array<std::array<double, 3>, 3> def_grad{};
    def_grad = this->ComputeDeltaDefGrad(nc, nenode, this->integrator_.alpha_f, dsf);

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
