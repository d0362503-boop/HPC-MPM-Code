#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::AssembleSolidSystem(const std::vector<double> &nvel_k,   //
                                              const std::vector<double> &naccel_k, //
                                              std::vector<std::array<double, 6>> &stress_k) {

    this->solid_.AssembleSystem(this->fsi_sys, naccel_k, nvel_k, stress_k);

    const double af = this->solid_.alpha_f;
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
