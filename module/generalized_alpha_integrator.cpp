#include <cmath>
#include <vector>

#include "module/generalized_alpha_integrator.h"

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
    this->dv_du = gamma / (beta * dt);
    this->dv_dv0 = gamma / beta - 1.0e0;
    this->dv_da0 = dt / 2.0e0 * (gamma / beta - 2.0e0);
    this->da_du = 1.0e0 / (beta * dt * dt);
    this->da_dv0 = 1.0e0 / (beta * dt);
    this->da_da0 = 1.0e0 / (2.0e0 * beta) - 1.0e0;

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
            this->dv_du * displacement[n] - this->dv_dv0 * old_velocity[n] - this->dv_da0 * old_acceleration[n];
        acceleration[n] =
            this->da_du * displacement[n] - this->da_dv0 * old_velocity[n] - this->da_da0 * old_acceleration[n];
    }

    return;
}
