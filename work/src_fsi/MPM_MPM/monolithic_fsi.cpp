#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <string>
#include <vector>

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::SolveFSISystem() {

    this->LumpedLagrangeMultiplier();

    this->BuildActiveDOFs();

    this->fluid_.MakeNSStabCoeff(this->fluid_.nvel);

    std::vector<std::array<double, 6>> stress_k = this->solid_.InitializeNRStress();

    std::vector<double> nvel_f(nodec * 3), nvel_s(nodec * 3);
    std::vector<double> naccel_f(nodec * 3), naccel_s(nodec * 3);
    std::array<double, 4> initial_norm{};

    VectorAssign(nodec * 3, this->fluid_.ndispl);
    VectorAssign(nodec * 3, this->solid_.ndispl);
    VectorAssign(nodec * 3, this->nlambda);
    VectorAssign(nodec, this->fluid_.npres);
    VectorAssign(nodec * 10, this->fsi_sys.x_lhs);

    for (int NR_it = 0; NR_it <= this->max_NR_it; NR_it++) {
        VectorAssign(this->fsi_sys.nmata, this->fsi_sys.amat);
        VectorAssign(nodec * 10, this->fsi_sys.b_rhs);

        this->fluid_.BCNRSet();
        this->fluid_.ComputeNodeVelAccelFromDispl(nvel_f, naccel_f);

        this->solid_.BCNRSet();
        this->solid_.ComputeNodeVelAccelFromDispl(nvel_s, naccel_s);

        this->AssembleFluidSystem(nvel_f, naccel_f);
        this->AssembleSolidSystem(nvel_s, naccel_s, stress_k);
        this->AssembleInterfaceSystem(nvel_f, nvel_s);

        int linear_iterations = this->SolveSystem(NR_it);

        this->UpdateNRIncrement();

        if (this->CheckNRConvergence(nvel_f, nvel_s, initial_norm, NR_it, linear_iterations)) { break; }
    }

    return;
}

void MPMMPMMonolithicFSI::AssembleInterfaceSystem(const std::vector<double> &nvel_f, //
                                                  const std::vector<double> &nvel_s) {

    for (int n = 0; n < nodec; n++) {
        int ncol = 0;
        int ida = this->fsi_sys.FindIndex(n, n, ncol);
        const double lm = this->nlm_lump_local[n];
        const double nb_para_f = this->fluid_.nb_para[0];
        const double nb_para_s = this->solid_.nb_para[0];
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[31]] -= nb_para_f * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[33]] -= nb_para_f * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[35]] -= nb_para_f * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[32]] += nb_para_s * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[34]] += nb_para_s * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[36]] += nb_para_s * lm;
    }

    for (int n = 0; n < nodec; n++) {
        this->fsi_sys.b_rhs[n + nodec * 7] = this->nlm_lump[n] * (nvel_f[n + nuc] - nvel_s[n + nuc]);
        this->fsi_sys.b_rhs[n + nodec * 8] = this->nlm_lump[n] * (nvel_f[n + nvc] - nvel_s[n + nvc]);
        this->fsi_sys.b_rhs[n + nodec * 9] = this->nlm_lump[n] * (nvel_f[n + nwc] - nvel_s[n + nwc]);
    }

    return;
}

void MPMMPMMonolithicFSI::AddLagrangeMultiplierToRHS(double af_coeff, const std::vector<int> &offsets) {

    for (int n = 0; n < nodec; n++) {
        this->fsi_sys.b_rhs[n + offsets[0]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + nuc];
        this->fsi_sys.b_rhs[n + offsets[1]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + nvc];
        this->fsi_sys.b_rhs[n + offsets[2]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + nwc];
    }

    return;
}

void MPMMPMMonolithicFSI::LumpedLagrangeMultiplier() {

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    double volp = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nlm_lump_local);
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
                    MakeSF(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

                    double phi_f = 0.0e0, phi_s = 0.0e0;
                    std::array<double, 3> grad_phi_f{}, grad_phi_s{};
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        double dsfi1 = dsf[ni][0];
                        double dsfi2 = dsf[ni][1];
                        double dsfi3 = dsf[ni][2];
                        phi_f += sfi * this->fluid_.nphi[nid];
                        phi_s += sfi * this->solid_.nphi[nid];
                        grad_phi_f[0] += dsfi1 * this->fluid_.nphi[nid];
                        grad_phi_f[1] += dsfi2 * this->fluid_.nphi[nid];
                        grad_phi_f[2] += dsfi3 * this->fluid_.nphi[nid];
                        grad_phi_s[0] += dsfi1 * this->solid_.nphi[nid];
                        grad_phi_s[1] += dsfi2 * this->solid_.nphi[nid];
                        grad_phi_s[2] += dsfi3 * this->solid_.nphi[nid];
                    }
                    std::array<double, 3> grad_phi;
                    for (int i = 0; i < 3; i++) { grad_phi[i] = phi_s * grad_phi_f[i] - phi_f * grad_phi_s[i]; }
                    double norm_grad_phi = NormVec3(grad_phi);

                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        this->nlm_lump_local[nid] += sfi * volp * norm_grad_phi;
                    }
                }
            }
        }
    }

    VectorAssign(nodec, this->nlm_lump);
    this->nlm_lump = this->nlm_lump_local;

    NodeVarComm(this->nlm_lump, 0);

    return;
}

void MPMMPMMonolithicFSI::UpdateNRIncrement() {

    for (int n = 0; n < nodec * 3; n++) {
        this->fluid_.ndispl[n] += this->fsi_sys.x_lhs[n];
        this->solid_.ndispl[n] += this->fsi_sys.x_lhs[n + nodec * 4];
        this->nlambda[n] += this->fsi_sys.x_lhs[n + nodec * 7];
    }

    for (int n = 0; n < nodec; n++) { this->fluid_.npres[n] += this->fsi_sys.x_lhs[n + npc]; }

    return;
}

bool MPMMPMMonolithicFSI::CheckNRConvergence(const std::vector<double> &nvel_f, const std::vector<double> &nvel_s,
                                             std::array<double, 4> &initial_norm, int NR_it, int linear_iterations) {

    const std::array<double, 4> absolute_tol = {1.0e-8, 1.0e-10, 1.0e-8, 1.0e-8};

    // stats[0]: fluid momentum residual [N].
    // stats[1]: continuity/PSPG residual [m^3/s].
    // stats[2]: solid momentum residual [N].
    // stats[3]: endpoint velocity difference [m/s].
    // stats[4..7]: corresponding active-component counts.
    // Accumulate squares below, then convert to global RMS; field RHS is unscaled.
    std::array<double, 8> stats{};

    // jump_max[0]: maximum velocity difference [m/s].
    // jump_max[1]: maximum displacement-increment difference [m].
    // Global componentwise maxima for diagnostics only; neither controls convergence.
    std::array<double, 2> jump_max{};

    for (int n : this->fsi_sys.owned_natural_ids) {
        for (int d = 0; d < 10; d++) {
            if (this->fixed_dof[n + d * nodec]) continue;
            const int field = d < 3 ? 0 : (d == 3 ? 1 : (d < 7 ? 2 : 3));
            const double residual =
                d < 7 ? this->fsi_sys.b_rhs[n + d * nodec] : nvel_s[n + (d - 7) * nodec] - nvel_f[n + (d - 7) * nodec];
            stats[field] += residual * residual;
            stats[field + 4] += 1.0;
            if (d >= 7) {
                jump_max[0] = std::max(jump_max[0], std::abs(residual));
                jump_max[1] = std::max(jump_max[1], std::abs(this->solid_.ndispl[n + (d - 7) * nodec] -
                                                             this->fluid_.ndispl[n + (d - 7) * nodec]));
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, stats.data(), 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, jump_max.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // All four RMS residuals must be finite and satisfy their thresholds.
    // Fields 0..2: max(absolute_tol, 1e-4 * this step's initial RMS).
    // Field 3: n+1 velocity-difference RMS <= 1e-8 m/s, without interface weights.
    bool converged = true;
    for (int field = 0; field < 4; field++) {
        stats[field] = std::sqrt(stats[field] / std::max(1.0, stats[field + 4]));
        if (NR_it == 0) { initial_norm[field] = stats[field]; }
        const double tolerance = field == 3 ? absolute_tol[field] : std::max(absolute_tol[field], 1.0e-4 * initial_norm[field]);
        converged = converged && stats[field] <= tolerance;
    }

    if (converged) {
        if (myrank == 0) {
            std::cout << "Monolithic_converge: " << std::setw(15) << NR_it << std::setw(15) << linear_iterations //
                      << std::scientific << std::setw(15) << stats[0] << std::setw(15) << stats[1]               //
                      << std::setw(15) << stats[2] << std::setw(15) << stats[3] << "\n";
        }
        return true;
    }

    if (myrank == 0) {
        std::cout << "Monolithic_NR: " << std::setw(15) << NR_it << std::scientific << std::setw(15) << stats[0] //
                  << std::setw(15) << stats[1] << std::setw(15) << stats[2] << std::setw(15) << stats[3] << "\n";
    }

    if (NR_it == this->max_NR_it) {
        if (myrank == 0) { std::cout << "Monolithic Newton did not converge" << "\n"; }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return false;
}
