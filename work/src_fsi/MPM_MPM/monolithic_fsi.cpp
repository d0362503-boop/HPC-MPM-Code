#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

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

    int solver_it = 0;
    for (int NR_it = 0; NR_it <= this->max_NR_it; NR_it++) {
        VectorAssign(this->fsi_sys.nmata, this->fsi_sys.amat);
        VectorAssign(nodec * 10, this->fsi_sys.b_rhs);

        this->fluid_.BCNRSet();
        this->fluid_.ComputeNodeVelAccelFromDispl(nvel_f, naccel_f);

        this->solid_.BCNRSet();
        this->solid_.ComputeNodeVelAccelFromDispl(nvel_s, naccel_s);

        this->AssembleFluidSystem(nvel_f, naccel_f);
        this->AssembleSolidSystem(nvel_s, naccel_s, stress_k);
        this->AssembleInterfaceSystem();

        if (this->CheckNRConvergence(initial_norm, NR_it, solver_it)) { break; }

        solver_it = this->SolveSystem(NR_it);

        this->UpdateNRIncrement();
    }

    return;
}

void MPMMPMMonolithicFSI::AssembleFluidSystem(const std::vector<double> &nvel_k, //
                                              const std::vector<double> &naccel_k) {

    this->fluid_.AssembleSystem(this->fsi_sys, nvel_k, naccel_k);

    const double af = this->fluid_.integrator_.alpha_f;
    for (int n = 0; n < nodec; n++) {
        int ncol = 0;
        int ida = this->fsi_sys.FindIndex(n, n, ncol);
        double lm = this->nlm_lump_local[n];
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[4]] -= af * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[9]] -= af * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[14]] -= af * lm;
    }

    this->AddLagrangeMultiplierToRHS(af, this->fluid_.RHSOffsets());

    return;
}

void MPMMPMMonolithicFSI::AssembleSolidSystem(const std::vector<double> &nvel_k,   //
                                              const std::vector<double> &naccel_k, //
                                              std::vector<std::array<double, 6>> &stress_k) {

    this->solid_.AssembleSystem(this->fsi_sys, naccel_k, nvel_k, stress_k);

    const double af = this->solid_.integrator_.alpha_f;
    for (int n = 0; n < nodec; n++) {
        int ncol = 0;
        int ida = this->fsi_sys.FindIndex(n, n, ncol);
        double lm = this->nlm_lump_local[n];
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[22]] += af * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[26]] += af * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[30]] += af * lm;
    }

    this->AddLagrangeMultiplierToRHS(-af, this->solid_.RHSOffsets());

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

std::array<double, 4> MPMMPMMonolithicFSI::ComputeNRResidualNorms() const {

    // stats[0]: fluid momentum residual [N].
    // stats[1]: continuity residual [m^3/s].
    // stats[2]: solid momentum residual [N].
    // stats[3]: displacement increment difference [m].
    // stats[4..7]: corresponding active-component counts.
    std::array<double, 8> stats{};

    for (int n : this->fsi_sys.owned_natural_ids) {
        for (int d = 0; d < 10; d++) {
            const int index = n + d * nodec;
            if (this->fixed_dof[index]) continue;

            int field;
            if (d < 3) { // Fluid momentum
                field = 0;
            } else if (d == 3) { // Fluid continuity
                field = 1;
            } else if (d < 7) { // Solid momentum
                field = 2;
            } else { // Interface displacement
                field = 3;
            }

            double residual = this->fsi_sys.b_rhs[index];
            if (field == 3) {
                const int index = n + (d - 7) * nodec;
                residual = this->solid_.ndispl[index] - this->fluid_.ndispl[index];
            }
            stats[field] += residual * residual;
            stats[field + 4] += 1.0;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, stats.data(), 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    std::array<double, 4> norms;
    for (int field = 0; field < 4; field++) { norms[field] = std::sqrt(stats[field] / std::max(1.0, stats[field + 4])); }

    return norms;
}

bool MPMMPMMonolithicFSI::CheckNRConvergence(std::array<double, 4> &initial_norm, int NR_it, int solver_it) {

    const std::array<double, 4> absolute_tol = {1.0e-8, 1.0e-10, 1.0e-8, 1.0e-10};
    const std::array<double, 4> stats = this->ComputeNRResidualNorms();

    // All four RMS residuals must be finite and satisfy their thresholds.
    // Fields 0..2: max(absolute_tol, 1e-4 * initial nonlinear residual RMS).
    // Field 3: displacement-increment RMS <= 1e-8 m, without interface weights.
    bool converged = true;
    for (int field = 0; field < 4; field++) {
        if (NR_it == 0 && field < 3) { initial_norm[field] = stats[field]; }
        const double tolerance = field == 3 ? absolute_tol[field] : std::max(absolute_tol[field], 1.0e-4 * initial_norm[field]);
        converged = converged && std::isfinite(stats[field]) && std::isfinite(initial_norm[field]) && stats[field] <= tolerance;
    }

    if (converged) {
        if (myrank == 0) {
            std::cout << "Monolithic_converge: " << std::setw(15) << NR_it << std::setw(15) << solver_it //
                      << std::scientific << std::setw(15) << stats[0] << std::setw(15) << stats[1]       //
                      << std::setw(15) << stats[2] << std::setw(15) << stats[3] << "\n";
        }
        return true;
    }

    if (myrank == 0) {
        std::cout << "Monolithic_NR: " << std::setw(15) << NR_it << std::setw(15) << solver_it //
                  << std::scientific << std::setw(15) << stats[0] << std::setw(15) << stats[1] //
                  << std::setw(15) << stats[2] << std::setw(15) << stats[3] << "\n";
    }

    if (NR_it == this->max_NR_it) {
        if (myrank == 0) { std::cout << "Monolithic Newton did not converge" << "\n"; }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return false;
}
