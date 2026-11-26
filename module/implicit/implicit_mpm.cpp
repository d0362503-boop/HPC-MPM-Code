#include "module/implicit/implicit_mpm.h"

#include <algorithm>
#include <vector>

#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"

void ImplicitMPM::ComputeNodeVelAccelFromDispl(const MaterialPoint &mp, std::vector<double> &nvel_k,
                                               std::vector<double> &naccel_k) const noexcept {

    this->integrator_.ComputeNodeVelAccelFromDispl(mp.ndispl, mp.nvel, mp.naccel, nvel_k, naccel_k);

    return;
}

double ImplicitMPM::LumpedMass(const MaterialPoint &mp, int pid, double sfi) const noexcept {
    double emd = this->integrator_.da_du * sfi * mp.mass[pid];

    return emd;
}

void ImplicitMPM::AddInertialForce(const MaterialPoint &mp, CrsMat &mat, const std::vector<double> &naccel, //
                                   const std::vector<int> &offsets) const {
    const ControlMesh &cmesh = mp.grid_.cmesh;

    for (int n = 0; n < cmesh.nodec; n++) {
        // --- For Generalized-α (if α_m = 1, back to Newmark-β) ---
        mat.b_rhs[n + offsets[0]] -= mp.nmass[n] * naccel[n + cmesh.nuc];
        mat.b_rhs[n + offsets[1]] -= mp.nmass[n] * naccel[n + cmesh.nvc];
        mat.b_rhs[n + offsets[2]] -= mp.nmass[n] * naccel[n + cmesh.nwc];
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
                mp.vel[n][i] += mp.params_.dt * ((1.0e0 - this->integrator_.gamma_nb) * accel_old[n][i] //
                                                 + this->integrator_.gamma_nb * mp.accel[n][i]);
            }
            mp.coord[n][i] += disp[n][i];
            if (has_shift) mp.coord[n][i] += disp_corr[n][i];
        }
    }

    return;
}
