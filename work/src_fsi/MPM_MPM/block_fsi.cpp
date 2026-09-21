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

#include "work/src_fsi/MPM_MPM/block_fsi.h"

using namespace mpmmpmblockfsi;

MPMMPMBlockFSI::~MPMMPMBlockFSI() {
    KSPDestroy(&this->schur_fluid_ksp_);
    KSPDestroy(&this->schur_solid_ksp_);
}

PetscErrorCode MPMMPMBlockFSI::ApplyExactSchurShell(Mat mat, Vec trial, Vec response) {
    void *context = nullptr;
    MatShellGetContext(mat, &context);
    return static_cast<MPMMPMBlockFSI *>(context)->ApplyExactSchur(trial, response);
}

PetscErrorCode MPMMPMBlockFSI::ApplyApproximateSchurShell(Mat mat, Vec trial, Vec response) {
    void *context = nullptr;
    MatShellGetContext(mat, &context);
    return static_cast<MPMMPMBlockFSI *>(context)->ApplyApproximateSchur(trial, response);
}

PetscErrorCode MPMMPMBlockFSI::SolveApproximateSchur(PC pc, Vec load, Vec solution) {
    void *context = nullptr;
    PCShellGetContext(pc, &context);
    return KSPSolve(static_cast<KSP>(context), load, solution);
}

void MPMMPMBlockFSI::BuildFluidResponse(int block_it) {
    MatDuplicate(this->fluid_.NS_.petsc_mat, MAT_COPY_VALUES, &this->schur_fluid_response_matrix_);
    MatEliminateZeros(this->schur_fluid_response_matrix_, PETSC_TRUE);

    if (block_it > 0) {
        KSPSetOperators(this->schur_fluid_ksp_, this->schur_fluid_response_matrix_, this->schur_fluid_response_matrix_);
        KSPSetReusePreconditioner(this->schur_fluid_ksp_, PETSC_TRUE);
        return;
    }

    KSPDestroy(&this->schur_fluid_ksp_);
    KSPCreate(PETSC_COMM_WORLD, &this->schur_fluid_ksp_);
    KSPSetOptionsPrefix(this->schur_fluid_ksp_, "fsi_response_");
    KSPSetType(this->schur_fluid_ksp_, KSPFGMRES);
    KSPGMRESSetRestart(this->schur_fluid_ksp_, 60);
    KSPSetTolerances(this->schur_fluid_ksp_, 1.0e-8, 1.0e-16, 1.0e6, 1000);
    KSPSetOperators(this->schur_fluid_ksp_, this->schur_fluid_response_matrix_, this->schur_fluid_response_matrix_);

    PC field_split_pc;
    KSPGetPC(this->schur_fluid_ksp_, &field_split_pc);
    PCSetType(field_split_pc, PCFIELDSPLIT);

    const PetscInt velocity_fields[] = {0, 1, 2};
    const PetscInt pressure_field = 3;
    PCFieldSplitSetBlockSize(field_split_pc, 4);
    PCFieldSplitSetFields(field_split_pc, "velocity", 3, velocity_fields, velocity_fields);
    PCFieldSplitSetFields(field_split_pc, "pressure", 1, &pressure_field, &pressure_field);
    PCFieldSplitSetType(field_split_pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(field_split_pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(field_split_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);
    KSPSetUp(this->schur_fluid_ksp_);

    PetscInt count;
    KSP *response_splits;
    PCFieldSplitGetSubKSP(field_split_pc, &count, &response_splits);
    for (int i = 0; i < count; i++) {
        PC response_pc;
        Mat response_matrix;

        KSPSetType(response_splits[i], KSPPREONLY);
        KSPGetPC(response_splits[i], &response_pc);
        KSPGetOperators(response_splits[i], nullptr, &response_matrix);
        MatEliminateZeros(response_matrix, PETSC_TRUE);
        PCSetType(response_pc, PCHYPRE);
        PCHYPRESetType(response_pc, "boomeramg");
        const std::string prefix = i == 0 ? "fsi_velocity_" : "fsi_pressure_";
        PCSetOptionsPrefix(response_pc, prefix.c_str());

        PetscOptionsSetValue(nullptr, ("-" + prefix + "pc_hypre_boomeramg_coarsen_type").c_str(), "hmis");
        PetscOptionsSetValue(nullptr, ("-" + prefix + "pc_hypre_boomeramg_interp_type").c_str(), "ext+i");
        PetscOptionsSetValue(nullptr, ("-" + prefix + "pc_hypre_boomeramg_relax_type_all").c_str(), "l1scaled-SOR/Jacobi");
        PetscOptionsSetValue(nullptr, ("-" + prefix + "pc_hypre_boomeramg_strong_threshold").c_str(), "0.7");
        if (i == 1) {
            PetscOptionsSetValue(nullptr, "-fsi_pressure_pc_hypre_boomeramg_smooth_type", "ILU");
            PetscOptionsSetValue(nullptr, "-fsi_pressure_pc_hypre_boomeramg_smooth_num_levels", "1");
            PetscOptionsSetValue(nullptr, "-fsi_pressure_pc_hypre_boomeramg_ilu_level", "0");
        }
        PCSetFromOptions(response_pc);
    }

    PetscFree(response_splits);
}

PetscErrorCode MPMMPMBlockFSI::ApplyApproximateSchur(Vec trial, Vec response) {
    PC fluid_pc, solid_pc;
    KSPGetPC(this->schur_fluid_ksp_, &fluid_pc);
    KSPGetPC(this->schur_solid_ksp_, &solid_pc);

    MatMult(this->schur_fluid_coupling_, trial, this->schur_fluid_rhs_);
    PCApply(fluid_pc, this->schur_fluid_rhs_, this->schur_fluid_solution_);
    for (int sweep = 0; sweep < 3; sweep++) {
        MatMult(this->schur_fluid_response_matrix_, this->schur_fluid_solution_, this->schur_fluid_defect_);
        VecAYPX(this->schur_fluid_defect_, -1.0, this->schur_fluid_rhs_);
        PCApply(fluid_pc, this->schur_fluid_defect_, this->schur_fluid_correction_);
        VecAXPY(this->schur_fluid_solution_, 1.0, this->schur_fluid_correction_);
    }
    MatMultTranspose(this->schur_fluid_coupling_, this->schur_fluid_solution_, response);
    VecScale(response, this->fluid_.nb_para[0]);

    MatMult(this->schur_solid_coupling_, trial, this->schur_solid_rhs_);
    PCApply(solid_pc, this->schur_solid_rhs_, this->schur_solid_solution_);
    MatMultTranspose(this->schur_solid_coupling_, this->schur_solid_solution_, this->schur_solid_response_);
    VecAXPY(response, this->solid_.nb_para[0], this->schur_solid_response_);

    return PETSC_SUCCESS;
}

void MPMMPMBlockFSI::BuildSolidResponse(int block_it) {

    MatDuplicate(this->solid_.SM_.petsc_mat, MAT_COPY_VALUES, &this->schur_solid_response_matrix_);
    MatEliminateZeros(this->schur_solid_response_matrix_, PETSC_TRUE);

    if (block_it > 0) {
        KSPSetOperators(this->schur_solid_ksp_, this->schur_solid_response_matrix_, this->schur_solid_response_matrix_);
        KSPSetReusePreconditioner(this->schur_solid_ksp_, PETSC_TRUE);
        return;
    }

    KSPDestroy(&this->schur_solid_ksp_);
    KSPCreate(PETSC_COMM_WORLD, &this->schur_solid_ksp_);
    KSPSetOperators(this->schur_solid_ksp_, this->schur_solid_response_matrix_, this->schur_solid_response_matrix_);
    KSPSetType(this->schur_solid_ksp_, KSPFGMRES);
    KSPSetOptionsPrefix(this->schur_solid_ksp_, "fsi_solid_response_");
    KSPSetTolerances(this->schur_solid_ksp_, 1.0e-8, 1.0e-16, 1.0e6, 1000);

    PC solid_pc;
    KSPGetPC(this->schur_solid_ksp_, &solid_pc);
    PCSetType(solid_pc, PCHYPRE);
    PCHYPRESetType(solid_pc, "boomeramg");
    PCSetOptionsPrefix(solid_pc, "fsi_solid_response_");
    PetscOptionsSetValue(nullptr, "-fsi_solid_response_pc_hypre_boomeramg_smooth_type", "Euclid");
    PetscOptionsSetValue(nullptr, "-fsi_solid_response_pc_hypre_boomeramg_smooth_num_levels", "1");
    PetscOptionsSetValue(nullptr, "-fsi_solid_response_pc_hypre_boomeramg_eu_level", "1");
    PCSetFromOptions(solid_pc);

    KSPSetUp(this->schur_solid_ksp_);

    return;
}

PetscInt MPMMPMBlockFSI::SolveResponse(Mat coupling, KSP solver, Vec trial_multiplier, Vec load, Vec solution, Vec response) {

    MatMult(coupling, trial_multiplier, load);

    KSPSetInitialGuessNonzero(solver, PETSC_FALSE);
    KSPSetReusePreconditioner(solver, PETSC_TRUE);
    KSPSolve(solver, load, solution);

    PetscInt iterations;
    KSPGetIterationNumber(solver, &iterations);
    KSPConvergedReason reason;
    KSPGetConvergedReason(solver, &reason);
    if (reason < 0 && myrank == 0) { std::cout << "PETSc KSP diverged, reason: " << reason << std::endl; }

    MatMultTranspose(coupling, solution, response);
    return iterations;
}

void MPMMPMBlockFSI::AddMultiplierIncrement(Vec delta_multiplier, const std::vector<PetscInt> &local_interface_ids) {
    std::vector<PetscInt> multiplier_rows(3 * this->fsi_intf.ibc);
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        for (int var = 0; var < 3; var++) { multiplier_rows[3 * n + var] = 3 * local_interface_ids[n] + var; }
    }

    IS from_is = nullptr, to_is = nullptr;
    Vec local_multiplier = nullptr;
    VecScatter multiplier_scatter = nullptr;
    ISCreateGeneral(PETSC_COMM_SELF, multiplier_rows.size(), multiplier_rows.data(), PETSC_COPY_VALUES, &from_is);
    ISCreateStride(PETSC_COMM_SELF, multiplier_rows.size(), 0, 1, &to_is);
    VecCreateSeq(PETSC_COMM_SELF, multiplier_rows.size(), &local_multiplier);
    VecScatterCreate(delta_multiplier, from_is, local_multiplier, to_is, &multiplier_scatter);
    VecScatterBegin(multiplier_scatter, delta_multiplier, local_multiplier, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(multiplier_scatter, delta_multiplier, local_multiplier, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *delta_values = nullptr;
    VecGetArrayRead(local_multiplier, &delta_values);
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        this->nfsi_force[nid + nuc] += PetscRealPart(delta_values[3 * n]);
        this->nfsi_force[nid + nvc] += PetscRealPart(delta_values[3 * n + 1]);
        this->nfsi_force[nid + nwc] += PetscRealPart(delta_values[3 * n + 2]);
    }
    VecRestoreArrayRead(local_multiplier, &delta_values);

    VecScatterDestroy(&multiplier_scatter);
    VecDestroy(&local_multiplier);
    ISDestroy(&from_is);
    ISDestroy(&to_is);
}

void MPMMPMBlockFSI::DestroySchurWorkspace() {
    MatDestroy(&this->schur_fluid_response_matrix_);
    MatDestroy(&this->schur_solid_response_matrix_);
    MatDestroy(&this->schur_fluid_coupling_);
    MatDestroy(&this->schur_solid_coupling_);

    VecDestroy(&this->schur_fluid_rhs_);
    VecDestroy(&this->schur_fluid_solution_);
    VecDestroy(&this->schur_fluid_defect_);
    VecDestroy(&this->schur_fluid_correction_);
    VecDestroy(&this->schur_solid_rhs_);
    VecDestroy(&this->schur_solid_solution_);
    VecDestroy(&this->schur_solid_response_);
}

void MPMMPMBlockFSI::BuildCouplingOperator(CrsMat &mat, const std::vector<double> &weights,
                                           const std::vector<PetscInt> &local_interface_ids, PetscInt local_multiplier_dofs,
                                           PetscInt global_multiplier_dofs, Mat &coupling, Vec &load, Vec &solution) {

    PetscInt local_field_dofs, global_field_dofs;
    VecGetLocalSize(mat.petsc_b, &local_field_dofs);
    VecGetSize(mat.petsc_b, &global_field_dofs);
    MatCreateAIJ(PETSC_COMM_WORLD, local_field_dofs, local_multiplier_dofs, global_field_dofs, global_multiplier_dofs, 1, nullptr,
                 1, nullptr, &coupling);

    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (mat.natural_is_owned[nid] == 0) continue;

        for (int var = 0; var < 3; var++) {
            const PetscInt row = mat.natural_var_gids[nid + var * nodec];
            if (std::binary_search(mat.petsc_bc_gids.begin(), mat.petsc_bc_gids.end(), row)) continue;
            const PetscInt column = 3 * local_interface_ids[n] + var;
            MatSetValue(coupling, row, column, weights[nid], INSERT_VALUES);
        }
    }
    MatAssemblyBegin(coupling, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(coupling, MAT_FINAL_ASSEMBLY);
    VecDuplicate(mat.petsc_b, &load);
    VecDuplicate(mat.petsc_x, &solution);

    return;
}

PetscErrorCode MPMMPMBlockFSI::ApplyExactSchur(Vec trial_multiplier, Vec response) {

    this->schur_fluid_iterations_ += this->SolveResponse(this->schur_fluid_coupling_, this->schur_fluid_ksp_, trial_multiplier,
                                                         this->schur_fluid_rhs_, this->schur_fluid_solution_, response);
    VecScale(response, this->fluid_.nb_para[0]);

    this->schur_solid_iterations_ +=
        this->SolveResponse(this->schur_solid_coupling_, this->schur_solid_ksp_, trial_multiplier, this->schur_solid_rhs_,
                            this->schur_solid_solution_, this->schur_solid_response_);
    VecAXPY(response, this->solid_.nb_para[0], this->schur_solid_response_);

    return PETSC_SUCCESS;
}

void MPMMPMBlockFSI::DetectFSIInterface() {

    VectorAssign(nodec, this->fsi_intf.nbc);

    this->LumpedLagrangeMultiplier();
    this->fluid_.NS_.BuildActiveRowMask();
    this->solid_.SM_.BuildActiveRowMask();

    this->fsi_intf.ibc = 0;
    for (int n = 0; n < nodec; n++) {
        if (this->fluid_.lm_lumped[n] > mtol && this->fluid_.NS_.active_row_mask[n] != 0 &&
            this->solid_.SM_.active_row_mask[n] != 0) {
            this->fsi_intf.nbc[this->fsi_intf.ibc++] = n;
        }
    }

    VectorAssign(nodec * 3, this->nfsi_force);

    return;
}

void MPMMPMBlockFSI::LumpedLagrangeMultiplier() {

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    double volp = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->fluid_.lm_lumped);
    VectorAssign(nodec, this->solid_.lm_lumped);
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
                        this->fluid_.lm_lumped[nid] += sfi * volp * norm_grad_phi;
                        this->solid_.lm_lumped[nid] += sfi * volp * norm_grad_phi;
                    }
                }
            }
        }
    }

    NodeVarComm(this->fluid_.lm_lumped, 0);
    NodeVarComm(this->solid_.lm_lumped, 0);

    return;
}

void MPMMPMBlockFSI::SolveFSISystem() {

    this->DetectFSIInterface();

    std::vector<double> fluid_velocity(nodec * 3), solid_velocity(nodec * 3), acceleration(nodec * 3);
    double rtr_fsi_ref, r0r_fsi_ref, rkr_fsi_ref, rtr_fsi_abs;
    for (int block_it = 0; block_it <= this->max_block_iter; block_it++) {

        this->fluid_.SolveNS();

        this->solid_.SolveSolid();

        this->fluid_.ComputeNodeVelAccelFromDispl(fluid_velocity, acceleration);

        this->solid_.ComputeNodeVelAccelFromDispl(solid_velocity, acceleration);

        if (block_it == 0) {
            this->CalFSIResidual(fluid_velocity, solid_velocity, r0r_fsi_ref, rtr_fsi_abs);
            rkr_fsi_ref = r0r_fsi_ref;
        } else {
            this->CalFSIResidual(fluid_velocity, solid_velocity, rkr_fsi_ref, rtr_fsi_abs);
        }
        rtr_fsi_ref = (r0r_fsi_ref > 1.0e-12) ? (rkr_fsi_ref / r0r_fsi_ref) : 0.0e0;

        // --- Convergence check ---
        if (rkr_fsi_ref < this->tol_ref * r0r_fsi_ref || //
            rtr_fsi_abs < this->tol_abs) {
            if (myrank == 0) {
                std::cout << "Block_converge:" << std::setw(10) << block_it  //
                          << std::setw(15) << std::scientific << rtr_fsi_ref //
                          << std::setw(15) << std::scientific << r0r_fsi_ref //
                          << std::setw(15) << std::scientific << rtr_fsi_abs << "\n";
            }
            break;
        } else if (block_it == this->max_block_iter) {
            if (myrank == 0) {
                std::cout << "Block_it:" << std::setw(10) << block_it << "\n";
                std::cout << "Block iteration can not converge" << "\n";
            }
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        this->UpdateFSIMultiplier(block_it, fluid_velocity, solid_velocity);
    }

    return;
}

void MPMMPMBlockFSI::CalFSIResidual(const std::vector<double> &fluid_velocity, const std::vector<double> &solid_velocity,
                                    double &rtr_ref, double &rtr_dof) {

    double norm = 0.0e0, intf_num = 0.0e0;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        int nid = this->fsi_intf.nbc[n];

        double dr1 = solid_velocity[nid + nuc] - fluid_velocity[nid + nuc];
        double dr2 = solid_velocity[nid + nvc] - fluid_velocity[nid + nvc];
        double dr3 = solid_velocity[nid + nwc] - fluid_velocity[nid + nwc];

        norm += dr1 * dr1 * dbc[nid + nuc]   //
                + dr2 * dr2 * dbc[nid + nvc] //
                + dr3 * dr3 * dbc[nid + nwc];

        intf_num += dbc[nid];
    }

    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &intf_num, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    rtr_ref = std::sqrt(norm);
    rtr_dof = std::sqrt(norm / intf_num);

    return;
}

void MPMMPMBlockFSI::BuildLumpedSchurPreconditioner(const std::vector<PetscInt> &interface_nodes, PetscInt local_dofs,
                                                    PetscInt global_dofs, Mat &preconditioner_mat) {

    MatCreateAIJ(PETSC_COMM_WORLD, local_dofs, local_dofs, global_dofs, global_dofs, 3, nullptr, 3, nullptr, &preconditioner_mat);
    MatSetBlockSize(preconditioner_mat, 3);
    MatSetOption(preconditioner_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    PetscInt row_start, row_end;
    MatGetOwnershipRange(this->fluid_.NS_.petsc_mat, &row_start, &row_end);
    for (PetscInt n = 0; n < static_cast<PetscInt>(interface_nodes.size()); n++) {
        const PetscInt node_id = interface_nodes[n];
        if (4 * node_id < row_start || 4 * node_id >= row_end) continue;

        std::array<PetscInt, 4> fluid_rows;
        std::array<PetscInt, 3> solid_rows;
        std::array<PetscScalar, 16> fluid_block;
        std::array<PetscScalar, 9> solid_block;
        for (int var = 0; var < 4; var++) { fluid_rows[var] = 4 * node_id + var; }
        for (int var = 0; var < 3; var++) { solid_rows[var] = 3 * node_id + var; }
        MatGetValues(this->fluid_.NS_.petsc_mat, 4, fluid_rows.data(), 4, fluid_rows.data(), fluid_block.data());
        MatGetValues(this->solid_.SM_.petsc_mat, 3, solid_rows.data(), 3, solid_rows.data(), solid_block.data());

        std::array<std::array<double, 3>, 3> fluid_velocity_block, solid_velocity_block;
        std::array<bool, 3> fluid_bc, solid_bc;
        double gf = 0.0, gs = 0.0;
        for (int i = 0; i < 3; i++) {
            const PetscInt column = 3 * n + i;
            PetscScalar fluid_weight, solid_weight;
            MatGetValues(this->schur_fluid_coupling_, 1, &fluid_rows[i], 1, &column, &fluid_weight);
            MatGetValues(this->schur_solid_coupling_, 1, &solid_rows[i], 1, &column, &solid_weight);
            fluid_bc[i] = (fluid_weight == 0.0);
            solid_bc[i] = (solid_weight == 0.0);
            gf = std::max(gf, PetscRealPart(fluid_weight));
            gs = std::max(gs, PetscRealPart(solid_weight));
            for (int j = 0; j < 3; j++) {
                fluid_velocity_block[i][j] =
                    PetscRealPart(fluid_block[4 * i + j] - fluid_block[4 * i + 3] * fluid_block[12 + j] / fluid_block[15]);
                solid_velocity_block[i][j] = PetscRealPart(solid_block[3 * i + j]);
            }
        }
        const auto fluid_inverse = InvMat3(fluid_velocity_block);
        const auto solid_inverse = InvMat3(solid_velocity_block);
        const double fluid_factor = this->fluid_.nb_para[0] * gf * gf;
        const double solid_factor = this->solid_.nb_para[0] * gs * gs;
        std::array<PetscScalar, 9> schur_block;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                schur_block[3 * i + j] = (fluid_bc[i] || fluid_bc[j] ? 0.0e0 : fluid_factor * fluid_inverse[i][j]) +
                                         (solid_bc[i] || solid_bc[j] ? 0.0e0 : solid_factor * solid_inverse[i][j]);
            }
            // Fully prescribed components have no multiplier response.
            if (fluid_bc[i] && solid_bc[i]) { schur_block[3 * i + i] = 1.0e0; }
        }
        const PetscInt row = 3 * n;
        const std::array<PetscInt, 3> rows{row, row + 1, row + 2};
        MatSetValues(preconditioner_mat, 3, rows.data(), 3, rows.data(), schur_block.data(), INSERT_VALUES);
    }

    MatAssemblyBegin(preconditioner_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(preconditioner_mat, MAT_FINAL_ASSEMBLY);

    return;
}

void MPMMPMBlockFSI::UpdateFSIMultiplier(int block_it, const std::vector<double> &fluid_velocity,
                                         const std::vector<double> &solid_velocity) {

    std::vector<PetscInt> local_global_nodes(this->fsi_intf.ibc);
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        local_global_nodes[n] = this->fluid_.NS_.natural_block_gids[this->fsi_intf.nbc[n]];
    }

    std::vector<int> recv_counts(nprocs), displacements(nprocs);
    const int local_count = static_cast<int>(local_global_nodes.size());
    MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, PETSC_COMM_WORLD);
    std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.begin() + 1);
    const int gathered_count = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);

    std::vector<PetscInt> interface_nodes(gathered_count);
    MPI_Allgatherv(local_global_nodes.data(), local_count, MPIU_INT, interface_nodes.data(), recv_counts.data(),
                   displacements.data(), MPIU_INT, PETSC_COMM_WORLD);
    std::sort(interface_nodes.begin(), interface_nodes.end());
    interface_nodes.erase(std::unique(interface_nodes.begin(), interface_nodes.end()), interface_nodes.end());

    const PetscInt interface_count = static_cast<PetscInt>(interface_nodes.size());
    const PetscInt global_dofs = 3 * interface_count;
    std::vector<PetscInt> local_interface_ids(this->fsi_intf.ibc);
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        local_interface_ids[n] =
            static_cast<PetscInt>(std::lower_bound(interface_nodes.begin(), interface_nodes.end(), local_global_nodes[n]) -
                                  interface_nodes.begin());
    }

    const PetscInt local_interface_count = interface_count / nprocs + (myrank < interface_count % nprocs ? 1 : 0);
    const PetscInt local_dofs = 3 * local_interface_count;
    Vec schur_rhs = nullptr, delta_multiplier = nullptr;
    VecCreateMPI(PETSC_COMM_WORLD, local_dofs, global_dofs, &schur_rhs);
    VecDuplicate(schur_rhs, &delta_multiplier);
    VecZeroEntries(schur_rhs);
    VecZeroEntries(delta_multiplier);

    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (this->fluid_.NS_.natural_is_owned[nid] == 0) continue;

        const PetscInt row = 3 * local_interface_ids[n];
        const std::array<PetscInt, 3> rows{row, row + 1, row + 2};
        std::array<PetscScalar, 3> values;
        const double gf = this->fluid_.lm_lumped[nid];
        const double gs = this->solid_.lm_lumped[nid];
        for (int var = 0; var < 3; var++) {
            const int id = nid + var * nodec;
            values[var] = gf * fluid_velocity[id] - gs * solid_velocity[id];
        }
        VecSetValues(schur_rhs, 3, rows.data(), values.data(), INSERT_VALUES);
    }
    VecAssemblyBegin(schur_rhs);
    VecAssemblyEnd(schur_rhs);

    this->schur_fluid_iterations_ = 0;
    this->schur_solid_iterations_ = 0;

    this->BuildCouplingOperator(this->fluid_.NS_, this->fluid_.lm_lumped, local_interface_ids, local_dofs, global_dofs,
                                this->schur_fluid_coupling_, this->schur_fluid_rhs_, this->schur_fluid_solution_);
    this->BuildCouplingOperator(this->solid_.SM_, this->solid_.lm_lumped, local_interface_ids, local_dofs, global_dofs,
                                this->schur_solid_coupling_, this->schur_solid_rhs_, this->schur_solid_solution_);
    VecDuplicate(schur_rhs, &this->schur_solid_response_);
    VecDuplicate(this->schur_fluid_rhs_, &this->schur_fluid_defect_);
    VecDuplicate(this->schur_fluid_rhs_, &this->schur_fluid_correction_);

    Mat preconditioner_mat = nullptr;
    this->BuildLumpedSchurPreconditioner(interface_nodes, local_dofs, global_dofs, preconditioner_mat);

    Mat schur_mat = nullptr;
    MatCreateShell(PETSC_COMM_WORLD, local_dofs, local_dofs, global_dofs, global_dofs, this, &schur_mat);
    MatShellSetOperation(schur_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(&MPMMPMBlockFSI::ApplyExactSchurShell));

    this->BuildSolidResponse(block_it);
    this->BuildFluidResponse(block_it);

    KSP schur_ksp = nullptr;
    KSPCreate(PETSC_COMM_WORLD, &schur_ksp);
    KSPSetOptionsPrefix(schur_ksp, "fsi_schur_");
    KSPSetOperators(schur_ksp, schur_mat, preconditioner_mat);
    KSPSetType(schur_ksp, KSPFGMRES);
    KSPGMRESSetRestart(schur_ksp, 50);
    KSPSetPCSide(schur_ksp, PC_RIGHT);
    KSPSetNormType(schur_ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(schur_ksp, 1.0e-8, 1.0e-20, 1.0e6, 200);

    PC schur_pc = nullptr;
    KSPGetPC(schur_ksp, &schur_pc);

    Mat approximate_mat = nullptr;
    MatCreateShell(PETSC_COMM_WORLD, local_dofs, local_dofs, global_dofs, global_dofs, this, &approximate_mat);
    MatShellSetOperation(approximate_mat, MATOP_MULT,
                         reinterpret_cast<void (*)(void)>(&MPMMPMBlockFSI::ApplyApproximateSchurShell));

    KSP approximate_ksp = nullptr;
    KSPCreate(PETSC_COMM_WORLD, &approximate_ksp);
    KSPSetOperators(approximate_ksp, approximate_mat, preconditioner_mat);
    KSPSetType(approximate_ksp, KSPFGMRES);
    KSPSetTolerances(approximate_ksp, 5.0e-3, 1.0e-30, 1.0e6, 30);

    PC approximate_pc;
    KSPGetPC(approximate_ksp, &approximate_pc);
    PCSetType(approximate_pc, PCPBJACOBI);

    PCSetType(schur_pc, PCSHELL);
    PCShellSetContext(schur_pc, approximate_ksp);
    PCShellSetApply(schur_pc, &MPMMPMBlockFSI::SolveApproximateSchur);
    KSPSetFromOptions(schur_ksp);

    KSPSolve(schur_ksp, schur_rhs, delta_multiplier);

    KSPConvergedReason reason;
    PetscReal residual_norm;
    KSPGetConvergedReason(schur_ksp, &reason);
    KSPGetResidualNorm(schur_ksp, &residual_norm);

    if (reason > 0) { this->AddMultiplierIncrement(delta_multiplier, local_interface_ids); }

    if (myrank == 0) {
        std::cout << "Schur_FGMRES_blockPC:" << std::setw(6) << block_it << std::setw(12) << this->schur_fluid_iterations_
                  << std::setw(12) << this->schur_solid_iterations_ << std::setw(15) << std::scientific << residual_norm << "\n";
    }

    KSPDestroy(&schur_ksp);
    KSPDestroy(&approximate_ksp);

    MatDestroy(&approximate_mat);
    MatDestroy(&preconditioner_mat);
    MatDestroy(&schur_mat);

    VecDestroy(&delta_multiplier);
    VecDestroy(&schur_rhs);

    this->DestroySchurWorkspace();

    return;
}
