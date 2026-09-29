#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"

void MaterialPoint::GeneralizedAlphaParaSet() {
    this->integrator_.GeneralizedAlphaParaSet();
}

std::vector<double> MaterialPoint::ComputeNodeAccelFromVel() const {
    return this->integrator_.ComputeNodeAccelFromVel(this->nvel, this->nvel_old, this->naccel, dt);
}

void MaterialPoint::NewmarkBetaParaSet() {
    this->integrator_.NewmarkBetaParaSet(dt);
}

void MaterialPoint::ComputeNodeVelAccelFromDispl(std::vector<double> &nvel_k, std::vector<double> &naccel_k) const noexcept {
    this->integrator_.ComputeNodeVelAccelFromDispl(this->ndispl, this->nvel, this->naccel, nvel_k, naccel_k);
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
