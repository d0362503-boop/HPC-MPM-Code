#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "module/dataset.h"
#include "module/generalized_alpha_integrator.h"
#include "module/material_point.h"
#include "module/mesh.h"

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

void MaterialPoint::CommitImplicitParticleKinematics(const std::vector<std::array<double, 3>> &accel_old,
                                                     const std::vector<std::array<double, 3>> &disp,
                                                     const std::vector<std::array<double, 3>> &disp_corr) {

    // const bool has_shift = !disp_corr.empty();
    for (int n = 0; n < this->num; n++) {
        std::array<double, 3> advected_coord;
        bool advected_outside = false;
        for (int i = 0; i < 3; i++) {
            if (this->solswitch == MapScheme::FLIP) {
                this->vel[n][i] += dt * ((1.0e0 - this->integrator_.gamma_nb) * accel_old[n][i] //
                                         + this->integrator_.gamma_nb * this->accel[n][i]);
            }
            advected_coord[i] = this->coord[n][i] + disp[n][i];
            this->coord[n][i] = advected_coord[i] + disp_corr[n][i];
            // if (advected_coord[i] < xyminw[i] || advected_coord[i] > xymaxw[i]) { advected_outside = true; }
        }

        // Reject shifting-induced boundary crossings.
        // if (has_shift && !advected_outside) {
        //     for (int i = 0; i < 3; i++) {
        //         double shifted_coord = advected_coord[i] + disp_corr[n][i];
        //         if (shifted_coord > xyminw[i] && shifted_coord < xymaxw[i]) { this->coord[n][i] = shifted_coord; }
        //     }
        // }
    }

    return;
}
