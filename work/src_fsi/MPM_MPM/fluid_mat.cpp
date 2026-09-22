#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::AssembleFluidSystem(const std::vector<double> &nvel_k, //
                                              const std::vector<double> &naccel_k) {

    this->fluid_.AssembleSystem(this->fsi_sys, nvel_k, naccel_k);

    const double af = this->fluid_.alpha_f;
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
