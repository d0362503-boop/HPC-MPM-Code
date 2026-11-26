#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::AssembleFluidSystem(const std::vector<double> &nvel_k, //
                                              const std::vector<double> &naccel_k, const TimeState &ts) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->fluid_.AssembleSystem(this->fsi_sys, nvel_k, naccel_k, ts);

    const double af = this->fluid_.implicit_.integrator_.alpha_f;
    for (int n = 0; n < cmesh.nodec; n++) {
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

void MPMMPMMonolithicFSI::AssembleSolidSystem(const std::vector<double> &naccel_k, //
                                              std::vector<std::array<double, 6>> &stress_k, const TimeState &ts) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->solid_.AssembleSystem(this->fsi_sys, naccel_k, stress_k, ts);

    const double af = this->solid_.implicit_.integrator_.alpha_f;
    for (int n = 0; n < cmesh.nodec; n++) {
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

void MPMMPMMonolithicFSI::AssembleInterfaceSystem(const std::vector<double> &nvel_f,
                                                 const std::vector<double> &nvel_s) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    const double cf = this->fluid_.implicit_.integrator_.dv_du;
    const double cs = this->solid_.implicit_.integrator_.dv_du;

    for (int n = 0; n < cmesh.nodec; n++) {
        int ncol = 0;
        int ida = this->fsi_sys.FindIndex(n, n, ncol);
        const double lm = this->nlm_lump_local[n];
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[31]] -= cf * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[33]] -= cf * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[35]] -= cf * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[32]] += cs * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[34]] += cs * lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[36]] += cs * lm;
    }

    for (int n = 0; n < cmesh.nodec; n++) {
        this->fsi_sys.b_rhs[n + cmesh.nodec * 7] = this->nlm_lump[n] * (nvel_f[n + cmesh.nuc] - nvel_s[n + cmesh.nuc]);
        this->fsi_sys.b_rhs[n + cmesh.nodec * 8] = this->nlm_lump[n] * (nvel_f[n + cmesh.nvc] - nvel_s[n + cmesh.nvc]);
        this->fsi_sys.b_rhs[n + cmesh.nodec * 9] = this->nlm_lump[n] * (nvel_f[n + cmesh.nwc] - nvel_s[n + cmesh.nwc]);
    }

    return;
}

void MPMMPMMonolithicFSI::AddLagrangeMultiplierToRHS(double af_coeff, const std::vector<int> &offsets) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    for (int n = 0; n < cmesh.nodec; n++) {
        this->fsi_sys.b_rhs[n + offsets[0]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + cmesh.nuc];
        this->fsi_sys.b_rhs[n + offsets[1]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + cmesh.nvc];
        this->fsi_sys.b_rhs[n + offsets[2]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + cmesh.nwc];
    }

    return;
}

void MPMMPMMonolithicFSI::AssemblePetscMat(CrsMat &mat) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    MatZeroEntries(mat.petsc_mat);
    mat.active_row_mask.assign(cmesh.nodec, 1);
    std::vector<PetscInt> columns;
    std::vector<PetscScalar> values;

    for (int n = 0; n < cmesh.nodec; n++) {
        for (int d = 0; d < mat.ndof; d++) {
            if (this->fixed_dof[n + d * cmesh.nodec]) continue;
            columns.clear();
            values.clear();
            for (int b = 0; b < mat.num_block; b++) {
                if (mat.block_row[b] != d) continue;
                for (int j = mat.matrow[n]; j < mat.matrow[n + 1]; j++) {
                    const double value = mat.amat[j + mat.block_id[b]];
                    if (value == 0.0) continue;
                    columns.push_back(mat.natural_var_gids[mat.matcolid[j] + mat.block_col[b] * cmesh.nodec]);
                    values.push_back(value);
                }
            }
            const PetscInt row = mat.natural_var_gids[n + d * cmesh.nodec];
            MatSetValues(mat.petsc_mat, 1, &row, columns.size(), columns.data(), values.data(), ADD_VALUES);
        }
    }

    MatAssemblyBegin(mat.petsc_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat.petsc_mat, MAT_FINAL_ASSEMBLY);
}
