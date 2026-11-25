#include "module/implicit_mpm.h"

#include <algorithm>
#include <vector>

#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"

void ImplicitMPM::EndpointKinematics(const MaterialPoint &mp, std::vector<double> &nvel_k,
                                     std::vector<double> &naccel_k) const noexcept {

    this->integrator_.ComputeNodeVelAccelFromDispl(mp.ndispl, mp.nvel, mp.naccel, nvel_k, naccel_k);

    return;
}

void ImplicitMPM::ApplyDisplacementBC(MaterialPoint &mp, std::vector<double> &ndispl) {
    mp.ubc.BCSetDispl(nuc, ndispl, this->integrator_, mp.nvel, mp.naccel);
    mp.vbc.BCSetDispl(nvc, ndispl, this->integrator_, mp.nvel, mp.naccel);
    mp.wbc.BCSetDispl(nwc, ndispl, this->integrator_, mp.nvel, mp.naccel);

    return;
}

double ImplicitMPM::LumpedMass(const MaterialPoint &mp, int pid, double sfi) const noexcept {
    double emd = this->integrator_.da_du * sfi * mp.mass[pid];

    return emd;
}

void ImplicitMPM::AddInertialForce(const MaterialPoint &mp, CrsMat &mat, const std::vector<double> &naccel, //
                                   const std::vector<int> &offsets) const {

    for (int n = 0; n < nodec; n++) {
        // --- For Generalized-α (if α_m = 1, back to Newmark-β) ---
        mat.b_rhs[n + offsets[0]] -= mp.nmass[n] * naccel[n + nuc];
        mat.b_rhs[n + offsets[1]] -= mp.nmass[n] * naccel[n + nvc];
        mat.b_rhs[n + offsets[2]] -= mp.nmass[n] * naccel[n + nwc];
    };

    return;
}

void ImplicitMPM::CommitNodalKinematics(MaterialPoint &mp, const std::vector<double> &nvel_k,
                                        const std::vector<double> &naccel_k) const {

    std::copy(nvel_k.begin(), nvel_k.end(), mp.nvel.begin());
    std::copy(naccel_k.begin(), naccel_k.end(), mp.naccel.begin());

    return;
}

void ImplicitMPM::CommitParticleKinematics(MaterialPoint &mp, const std::vector<std::array<double, 3>> &accel_old,
                                           const std::vector<std::array<double, 3>> &disp,
                                           const std::vector<std::array<double, 3>> &disp_corr) const {

    const bool has_shift = !disp_corr.empty();
    for (int n = 0; n < mp.num; n++) {
        for (int i = 0; i < 3; i++) {
            if (mp.solswitch == MapScheme::FLIP) {
                mp.vel[n][i] += dt * ((1.0e0 - this->integrator_.gamma_nb) * accel_old[n][i] //
                                       + this->integrator_.gamma_nb * mp.accel[n][i]);
            }
            mp.coord[n][i] += disp[n][i];
            if (has_shift) mp.coord[n][i] += disp_corr[n][i];
        }
    }

    return;
}
