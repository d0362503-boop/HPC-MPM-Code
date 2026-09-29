#include "module/generalized_alpha_integrator.h"

#include <cmath>

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
}

std::vector<double> GeneralizedAlphaIntegrator::ComputeNodeAccelFromVel(
    const std::vector<double> &velocity, const std::vector<double> &old_velocity,
    const std::vector<double> &old_acceleration, double dt) const {
    const double para1 = 1.0e0 / (this->gamma_nb * dt);
    const double para2 = (1.0e0 - this->gamma_nb) / this->gamma_nb;
    std::vector<double> acceleration(velocity.size());
    for (std::size_t n = 0; n < velocity.size(); ++n) {
        acceleration[n] = para1 * (velocity[n] - old_velocity[n]) - para2 * old_acceleration[n];
    }
    return acceleration;
}

void GeneralizedAlphaIntegrator::ComputeNodeVelAccelFromDispl(
    const std::vector<double> &displacement, const std::vector<double> &old_velocity,
    const std::vector<double> &old_acceleration, std::vector<double> &velocity,
    std::vector<double> &acceleration) const noexcept {
    for (std::size_t n = 0; n < displacement.size(); ++n) {
        velocity[n] = this->nb_para[0] * displacement[n] - this->nb_para[1] * old_velocity[n]
                      - this->nb_para[2] * old_acceleration[n];
        acceleration[n] = this->nb_para[3] * displacement[n] - this->nb_para[4] * old_velocity[n]
                          - this->nb_para[5] * old_acceleration[n];
    }
}
