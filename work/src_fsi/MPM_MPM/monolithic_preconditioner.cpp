#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::ConfigurePreconditioner(CrsMat &mat, PC pc) {
    KSPSetOptionsPrefix(mat.ksp, "fsi_");
    PCSetOptionsPrefix(pc, "fsi_");
    KSPGMRESSetRestart(mat.ksp, 200);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetBlockSize(pc, 1);
    if (mat.HasReducedLayout()) {
        PetscCallAbort(MPI_COMM_WORLD, mat.SetFieldSplitComponents(pc, "fluid", {0, 1, 2, 3}));
        PetscCallAbort(MPI_COMM_WORLD, mat.SetFieldSplitComponents(pc, "coupled", {4, 5, 6, 7, 8, 9}));
    }
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);

    // Fluid velocity-pressure preconditioner.
    const char *options[][2] = {{"-fsi_ksp_gmres_cgs_refinement_type", "refine_always"},
                                {"-fsi_fieldsplit_fluid_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_fluid_ksp_rtol", "0.5"},
                                {"-fsi_fieldsplit_fluid_ksp_max_it", "100"},
                                {"-fsi_fieldsplit_fluid_pc_type", "fieldsplit"},
                                {"-fsi_fieldsplit_fluid_pc_fieldsplit_type", "schur"},
                                {"-fsi_fieldsplit_fluid_pc_fieldsplit_schur_fact_type", "lower"},
                                {"-fsi_fieldsplit_fluid_pc_fieldsplit_schur_precondition", "selfp"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_pc_type", "asm"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_shift_amount", "1e-3"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_pc_type", "asm"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_levels", "1"},
                                // Reject tiny pressure ILU pivots.
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_zeropivot", "1e-8"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_shift_amount", "1e-3"},

                                // Joint solid-multiplier iteration.
                                {"-fsi_fieldsplit_coupled_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_pc_type", "ksp"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_rtol", "0.5"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_max_it", "100"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_gmres_restart", "100"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_type", "fieldsplit"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_fieldsplit_type", "schur"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_fieldsplit_schur_fact_type", "lower"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_fieldsplit_schur_precondition", "selfp"},

                                // Solid ASM with local ILU(1).
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_rtol", "0.1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_max_it", "100"},
                                // {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_gmres_restart", "200"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_pc_type", "asm"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_shift_amount", "1e-3"},

                                // Multiplier Schur preconditioner.
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_pc_type", "asm"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_factor_shift_amount", "1e-3"}};

    for (const auto &option : options) {
        PetscBool supplied;
        PetscOptionsHasName(nullptr, nullptr, option[0], &supplied);
        if (!supplied) { PetscOptionsSetValue(nullptr, option[0], option[1]); }
    }
}

PetscErrorCode MPMMPMMonolithicFSI::ConfigureReducedPreconditioner(KSP ksp) {
    PetscFunctionBeginUser;
    PC pc, fluid_pc, coupled_pc, inner_pc;
    KSP *outer, inner;
    PetscInt count;
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCFieldSplitGetSubKSP(pc, &count, &outer));

    PetscCall(KSPGetPC(outer[0], &fluid_pc));
    PetscCall(this->fsi_sys.SetFieldSplitComponents(fluid_pc, "0", {0, 1, 2}, {0, 1, 2, 3}));
    PetscCall(this->fsi_sys.SetFieldSplitComponents(fluid_pc, "1", {3}, {0, 1, 2, 3}));

    // Iterate on PETSc's assembled fluid-eliminated approximation.
    PetscCall(KSPGetPC(outer[1], &coupled_pc));
    PetscCall(PCSetUseAmat(coupled_pc, PETSC_FALSE));
    PetscCall(PCKSPGetKSP(coupled_pc, &inner));
    Mat coupled_preconditioner;
    PetscCall(KSPGetOperators(outer[1], nullptr, &coupled_preconditioner));
    PetscCall(KSPSetOperators(inner, coupled_preconditioner, coupled_preconditioner));
    PetscCall(KSPSetFromOptions(inner));
    PetscCall(KSPGetPC(inner, &inner_pc));
    PetscCall(this->fsi_sys.SetFieldSplitComponents(inner_pc, "solid", {4, 5, 6}, {4, 5, 6, 7, 8, 9}));
    PetscCall(this->fsi_sys.SetFieldSplitComponents(inner_pc, "lambda", {7, 8, 9}, {4, 5, 6, 7, 8, 9}));
    PetscCall(KSPSetUp(inner));
    PetscCall(PetscFree(outer));
    PetscFunctionReturn(PETSC_SUCCESS);
}
