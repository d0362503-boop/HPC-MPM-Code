#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dataset.h"
#include "material_point.h"
#include "mesh.h"

void MaterialPoint::GeneralizedAlphaParaSet() {
    double temp = 1.0e0 + this->spec_rad;
    this->alpha_f = 1.0e0 / temp;

    if (this->ode_order == 2) {
        // --- for second order PDE (Solid) ---
        this->alpha_m = (2.0e0 - this->spec_rad) / temp;
    } else if (this->ode_order == 1) {
        // --- for first order PDE (Fluid) ---
        this->alpha_m = 0.5e0 * ((3.0e0 - this->spec_rad) / temp);
    }

    temp = 1.0e0 - this->alpha_f + this->alpha_m;
    this->gamma_nb = temp - 0.5e0;
    this->beta_nb = 0.25e0 * pow(temp, 2);

    return;
}

std::vector<double> MaterialPoint::GeneralizedAlphaNodeAccelUpdate() {
    double para1 = 1.0e0 / (this->gamma_nb * dt);
    double para2 = (1.0e0 - this->gamma_nb) / this->gamma_nb;

    std::vector<double> naccel_k(nodec * 3);
    for (int n = 0; n < nodec * 3; n++) {
        naccel_k[n] = para1 * (this->nvel[n] - this->nvel_old[n]) //
                      - para2 * this->naccel[n];
    }

    return naccel_k;
}

void MaterialPoint::NewmarkBetaParaSet() {
    VectorAssign(6, this->nb_para);

    double gamma = this->gamma_nb;
    double beta = this->beta_nb;

    // ---- Newmark beta parameter ----
    this->nb_para[0] = gamma / (beta * dt);
    this->nb_para[1] = gamma / beta - 1.0e0;
    this->nb_para[2] = dt / 2.0e0 * (gamma / beta - 2.0e0);
    this->nb_para[3] = 1.0e0 / (beta * dt * dt);
    this->nb_para[4] = 1.0e0 / (beta * dt);
    this->nb_para[5] = 1.0e0 / (2.0e0 * beta) - 1.0e0;

    return;
}

void MaterialPoint::PredictNewmarkBetaVelAndAccel(std::vector<double> &nvel_k, std::vector<double> &naccel_k) {
    for (int n = 0; n < nodec * 3; n++) {
        nvel_k[n] = this->nb_para[0] * this->ndispl[n] //
                    - this->nb_para[1] * this->nvel[n] //
                    - this->nb_para[2] * this->naccel[n];
        naccel_k[n] = this->nb_para[3] * this->ndispl[n] //
                      - this->nb_para[4] * this->nvel[n] //
                      - this->nb_para[5] * this->naccel[n];
    }

    return;
}

void MaterialPoint::CommitNodalKinematics(const std::vector<double> &nvel_k, const std::vector<double> &naccel_k) {
    std::copy(nvel_k.begin(), nvel_k.end(), this->nvel.begin());
    std::copy(naccel_k.begin(), naccel_k.end(), this->naccel.begin());

    return;
}

void MaterialPoint::CommitParticleKinematics(const std::vector<std::array<double, 3>> &accel_old,
                                             const std::vector<std::array<double, 3>> &disp,
                                             const std::vector<std::array<double, 3>> &delta_corr) {
    const bool has_shift = !delta_corr.empty();
    for (int n = 0; n < this->num; n++) {
        for (int i = 0; i < 3; i++) {
            if (this->solswitch == MapScheme::FLIP) {
                this->vel[n][i] += dt * ((1.0e0 - this->gamma_nb) * accel_old[n][i] //
                                         + this->gamma_nb * this->accel[n][i]);
            }
            this->coord[n][i] += disp[n][i];
            // ---- particle shifting ----
            if (has_shift) this->coord[n][i] += delta_corr[n][i];
        }
    }

    return;
}