#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "module/dataset.h"
#include "module/generalized_alpha_integrator.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

void GeneralizedAlphaIntegrator::GeneralizedAlphaParaSet() {

    double temp = 1.0e0 + this->spec_rad;
    this->alpha_f = 1.0e0 / temp;

    if (this->ode_order == 2) {
        this->alpha_m = (2.0e0 - this->spec_rad) / temp;
    } else if (this->ode_order == 1) {
        this->alpha_m = 0.5e0 * ((3.0e0 - this->spec_rad) / temp);
    }

    temp = 1.0e0 - this->alpha_f + this->alpha_m;
    this->gamma_nb = temp - 0.5e0;
    this->beta_nb = 0.25e0 * std::pow(temp, 2);

    return;
}

void GeneralizedAlphaIntegrator::NewmarkBetaParaSet(double dt) {

    const double gamma = this->gamma_nb;
    const double beta = this->beta_nb;
    this->nb_para[0] = gamma / (beta * dt);
    this->nb_para[1] = gamma / beta - 1.0e0;
    this->nb_para[2] = dt / 2.0e0 * (gamma / beta - 2.0e0);
    this->nb_para[3] = 1.0e0 / (beta * dt * dt);
    this->nb_para[4] = 1.0e0 / (beta * dt);
    this->nb_para[5] = 1.0e0 / (2.0e0 * beta) - 1.0e0;

    return;
}

std::vector<double> GeneralizedAlphaIntegrator::ComputeNodeAccelFromVel(const std::vector<double> &velocity,
                                                                        const std::vector<double> &old_velocity,
                                                                        const std::vector<double> &old_acceleration,
                                                                        double dt) const {

    const double para1 = 1.0e0 / (this->gamma_nb * dt);
    const double para2 = (1.0e0 - this->gamma_nb) / this->gamma_nb;
    std::vector<double> acceleration(velocity.size());
    for (std::size_t n = 0; n < velocity.size(); ++n) {
        acceleration[n] = para1 * (velocity[n] - old_velocity[n]) - para2 * old_acceleration[n];
    }

    return acceleration;
}

void GeneralizedAlphaIntegrator::ComputeNodeVelAccelFromDispl(const std::vector<double> &displacement,
                                                              const std::vector<double> &old_velocity,
                                                              const std::vector<double> &old_acceleration,
                                                              std::vector<double> &velocity,
                                                              std::vector<double> &acceleration) const noexcept {

    for (std::size_t n = 0; n < displacement.size(); ++n) {
        velocity[n] =
            this->nb_para[0] * displacement[n] - this->nb_para[1] * old_velocity[n] - this->nb_para[2] * old_acceleration[n];
        acceleration[n] =
            this->nb_para[3] * displacement[n] - this->nb_para[4] * old_velocity[n] - this->nb_para[5] * old_acceleration[n];
    }

    return;
}

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

    const auto closed_face = this->DetermineClosedBoundaryFaces();
    const bool has_shift = !disp_corr.empty();
    for (int n = 0; n < this->num; n++) {
        for (int i = 0; i < 3; i++) {
            if (this->solswitch == MapScheme::FLIP) {
                this->vel[n][i] += dt * ((1.0e0 - this->integrator_.gamma_nb) * accel_old[n][i] //
                                         + this->integrator_.gamma_nb * this->accel[n][i]);
            }
            this->coord[n][i] += disp[n][i];
            if (has_shift) this->coord[n][i] += disp_corr[n][i];
            this->coord[n][i] = this->ConstrainClosedBoundaryCoordinate(this->coord[n][i], i, closed_face);
        }
    }

    return;
}
