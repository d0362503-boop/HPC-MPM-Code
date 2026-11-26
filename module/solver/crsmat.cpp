#include "module/solver/crsmat.h"

#include <petsc.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

void CrsMat::BuildCrsMat(int num_block) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->num_block = num_block;

    if (this->use_petsc) { this->ResetPetscSolver(); }

    this->nmat = 0;

    int nxp = cmesh.xynodec[0];
    int nyp = cmesh.xynodec[1];
    int nzp = cmesh.xynodec[2];

    int idimp2 = *std::max_element(cmesh.idimc.begin(), cmesh.idimc.end());

    int nmatptmp = cmesh.nodec * std::pow((idimp2 + idimp2 + 1), 3);
    std::vector<int> matcolidp(nmatptmp, -1);
    VectorAssign(cmesh.nodec + 1, this->matrow, -1);

    for (int icnzp = 0; icnzp < nzp; ++icnzp) {
        int isnzp = icnzp - cmesh.idimc[2];
        if (isnzp < 0) isnzp = 0;
        int ienzp = icnzp + cmesh.idimc[2];
        if (ienzp >= nzp) ienzp = nzp - 1;

        for (int icnyp = 0; icnyp < nyp; ++icnyp) {
            int isnyp = icnyp - cmesh.idimc[1];
            if (isnyp < 0) isnyp = 0;
            int ienyp = icnyp + cmesh.idimc[1];
            if (ienyp >= nyp) ienyp = nyp - 1;

            for (int icnxp = 0; icnxp < nxp; ++icnxp) {
                int isnxp = icnxp - cmesh.idimc[0];
                if (isnxp < 0) isnxp = 0;
                int ienxp = icnxp + cmesh.idimc[0];
                if (ienxp >= nxp) ienxp = nxp - 1;

                int idn = icnxp + nxp * icnyp + nxp * nyp * icnzp;
                this->matrow[idn] = this->nmat;

                for (int inz = isnzp; inz <= ienzp; ++inz) {
                    for (int iny = isnyp; iny <= ienyp; ++iny) {
                        for (int inx = isnxp; inx <= ienxp; ++inx) {
                            int idnmat = inx + nxp * iny + nxp * nyp * inz;
                            matcolidp[this->nmat] = idnmat;
                            this->nmat++;
                        }
                    }
                }
            }
        }
    }

    VectorAssign(this->nmat, this->matcolid);
    for (int i = 0; i < this->nmat; i++) { this->matcolid[i] = matcolidp[i]; }

    this->matrow[cmesh.nodec] = this->nmat;

    this->nmata = this->nmat * num_block;
    VectorAssign(num_block, this->block_id);
    for (int i = 0; i < num_block; i++) { this->block_id[i] = this->nmat * i; }
    VectorAssign(this->nmata, this->amat);

    if (this->use_petsc) { this->InitPetscSolver(); }

    return;
}

void CrsMat::ResetPetscSolver() {

    this->ResetReducedSystem();

    if (this->scatter_to_all) { VecScatterDestroy(&this->scatter_to_all); }
    if (this->seq_x) { VecDestroy(&this->seq_x); }
    if (this->ksp) { KSPDestroy(&this->ksp); }
    if (this->petsc_mat) { MatDestroy(&this->petsc_mat); }
    if (this->petsc_b) { VecDestroy(&this->petsc_b); }
    if (this->petsc_x) { VecDestroy(&this->petsc_x); }
    if (this->petsc_var_lgmap_) { ISLocalToGlobalMappingDestroy(&this->petsc_var_lgmap_); }
    if (this->petsc_block_lgmap_) { ISLocalToGlobalMappingDestroy(&this->petsc_block_lgmap_); }

    this->prev_ksp_its_ = -1;
    this->force_rebuild_next_ = false;

    return;
}

void CrsMat::BuildNaturalGlobalMaps(const std::vector<int> &node_l2g) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    VectorAssign(cmesh.nodec, this->natural_block_gids);
    VectorAssign(cmesh.nodec * this->ndof, this->natural_var_gids);

    for (int n = 0; n < cmesh.nodec; ++n) {
        const PetscInt global_node = static_cast<PetscInt>(node_l2g[n]);
        this->natural_block_gids[n] = global_node;
        for (int var = 0; var < this->ndof; ++var) { this->natural_var_gids[n + var * cmesh.nodec] = global_node * this->ndof + var; }
    }

    return;
}

void CrsMat::BuildOwnershipMetadata(const std::vector<int> &node_l2g) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->owned_natural_ids.clear();
    this->ghost_natural_ids.clear();

    std::vector<int> all_aelemmin(nprocs * 3);
    std::vector<int> all_aelemmax(nprocs * 3);
    MPI_Allgather(mesh.aelemmin.data(), 3, MPI_INT, all_aelemmin.data(), 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(mesh.aelemmax.data(), 3, MPI_INT, all_aelemmax.data(), 3, MPI_INT, MPI_COMM_WORLD);

    for (int n = 0; n < cmesh.nodec; ++n) {
        const int global_node = node_l2g[n];
        const int kg = global_node / (cmesh.xynodecw[0] * cmesh.xynodecw[1]);
        const int rem = global_node % (cmesh.xynodecw[0] * cmesh.xynodecw[1]);
        const int jg = rem / cmesh.xynodecw[0];
        const int ig = rem % cmesh.xynodecw[0];

        bool is_owned = true;
        for (int p = 0; p < nprocs; ++p) {
            if (p == myrank) continue;

            bool in_p = true;
            for (int d = 0; d < 3; ++d) {
                const int p_min = all_aelemmin[p * 3 + d];
                const int p_max = all_aelemmax[p * 3 + d] + cmesh.idimc[d];
                const int coord = (d == 0) ? ig : (d == 1) ? jg : kg;
                if (coord < p_min || coord > p_max) {
                    in_p = false;
                    break;
                }
            }

            if (in_p) {
                for (int d = 0; d < 3; ++d) {
                    const int my_min = mesh.aelemmin[d];
                    const int p_min = all_aelemmin[p * 3 + d];
                    if (p_min < my_min) {
                        is_owned = false;
                        break;
                    } else if (p_min > my_min) {
                        break;
                    }
                }
                if (!is_owned) break;
            }
        }

        if (is_owned) {
            this->owned_natural_ids.push_back(n);
        } else {
            this->ghost_natural_ids.push_back(n);
        }
    }

    this->local_node = static_cast<int>(this->owned_natural_ids.size());
    this->ghost_node = static_cast<int>(this->ghost_natural_ids.size());

    return;
}

void CrsMat::BuildPetscLocalMaps() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->petsc_local_to_natural.clear();
    this->petsc_local_to_natural.reserve(cmesh.nodec);

    for (int nid : this->owned_natural_ids) this->petsc_local_to_natural.push_back(nid);
    for (int nid : this->ghost_natural_ids) this->petsc_local_to_natural.push_back(nid);

    this->natural_to_petsc_local.assign(cmesh.nodec, -1);
    for (int pos = 0; pos < cmesh.nodec; ++pos) { this->natural_to_petsc_local[this->petsc_local_to_natural[pos]] = pos; }

    VectorAssign(cmesh.nodec, this->petsc_local_block_gids);
    VectorAssign(cmesh.nodec * this->ndof, this->petsc_local_var_gids);
    for (int pos = 0; pos < cmesh.nodec; ++pos) {
        const int natural_id = this->petsc_local_to_natural[pos];
        const PetscInt global_node = this->natural_block_gids[natural_id];
        this->petsc_local_block_gids[pos] = global_node;
        for (int var = 0; var < this->ndof; ++var) {
            this->petsc_local_var_gids[pos + var * cmesh.nodec] = global_node * this->ndof + var;
        }
    }

    return;
}

void CrsMat::CheckOwnershipMetadata() const {
    const ControlMesh &cmesh = this->grid_.cmesh;

    if (static_cast<int>(this->owned_natural_ids.size()) != this->local_node) MPI_Abort(MPI_COMM_WORLD, 1);
    if (static_cast<int>(this->ghost_natural_ids.size()) != this->ghost_node) MPI_Abort(MPI_COMM_WORLD, 1);
    if (static_cast<int>(this->petsc_local_to_natural.size()) != cmesh.nodec) MPI_Abort(MPI_COMM_WORLD, 1);

    for (int n = 0; n < cmesh.nodec; ++n) {
        if (this->natural_to_petsc_local[n] < 0) MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int global_nodec = cmesh.xynodecw[0] * cmesh.xynodecw[1] * cmesh.xynodecw[2];
    int total_local_node = 0;
    MPI_Allreduce(&this->local_node, &total_local_node, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (total_local_node != global_nodec) MPI_Abort(MPI_COMM_WORLD, 1);

    return;
}

void CrsMat::BuildLGMAP() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    std::vector<int> node_l2g(cmesh.nodec);

    for (int kl = 0; kl < cmesh.xynodec[2]; kl++) {
        for (int jl = 0; jl < cmesh.xynodec[1]; jl++) {
            for (int il = 0; il < cmesh.xynodec[0]; il++) {
                int kg = kl + mesh.aelemmin[2];
                int jg = jl + mesh.aelemmin[1];
                int ig = il + mesh.aelemmin[0];
                int nid_l = il + cmesh.xynodec[0] * jl + cmesh.xynodec[0] * cmesh.xynodec[1] * kl;
                int nid_g = ig + cmesh.xynodecw[0] * jg + cmesh.xynodecw[0] * cmesh.xynodecw[1] * kg;
                node_l2g[nid_l] = nid_g;
            }
        }
    }

    this->BuildNaturalGlobalMaps(node_l2g);
    BuildOwnershipMetadata(node_l2g);
    this->BuildPetscLocalMaps();
    CheckOwnershipMetadata();

    return;
}

void CrsMat::BuildPetscMat() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    PetscInt local_n = static_cast<PetscInt>(this->local_node) * this->ndof;

    MatCreate(MPI_COMM_WORLD, &this->petsc_mat);
    MatSetSizes(this->petsc_mat, local_n, local_n, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(this->petsc_mat, MATAIJ);
    MatSetBlockSize(this->petsc_mat, this->ndof);

    // Compute the maximum nonzeros per row across ALL nodec (owned + ghost).
    // Ghost rows are used during assembly (remote contributions), so the
    // reusable column/value buffers must be large enough for them too.
    int max_row_nnz = 0;
    for (int ii = 0; ii < cmesh.nodec; ++ii) {
        int row_nnz = this->matrow[ii + 1] - this->matrow[ii];
        if (row_nnz > max_row_nnz) max_row_nnz = row_nnz;
    }

    PetscInt prealloc = static_cast<PetscInt>(max_row_nnz) * this->ndof;
    MatMPIAIJSetPreallocation(this->petsc_mat, prealloc, NULL, prealloc, NULL);
    MatSeqAIJSetPreallocation(this->petsc_mat, prealloc, NULL);

    // Keep the sparsity pattern across solves so that repeated assemblies
    // do not reallocate internal PETSc data structures.
    MatSetOption(this->petsc_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    // Allow dynamic allocation for new nonzeros caused by moving material
    // points (MPM) or changing active sets, preventing PETSc malloc errors.
    MatSetOption(this->petsc_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    // Preserve diagonal storage slots.
    PetscInt row_begin, row_end;
    MatGetOwnershipRange(this->petsc_mat, &row_begin, &row_end);
    for (PetscInt row = row_begin; row < row_end; row++) { MatSetValue(this->petsc_mat, row, row, 0.0, INSERT_VALUES); }
    MatAssemblyBegin(this->petsc_mat, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(this->petsc_mat, MAT_FLUSH_ASSEMBLY);

    // Create distributed RHS and solution vectors matching the matrix layout.
    VecCreateMPI(MPI_COMM_WORLD, local_n, PETSC_DETERMINE, &this->petsc_b);
    VecCreateMPI(MPI_COMM_WORLD, local_n, PETSC_DETERMINE, &this->petsc_x);

    // Bind the local-to-global mapping so that all local indices used in
    // MatSetValuesBlockedLocal / VecSetValuesLocal are interpreted with the
    // owned-first PETSc-local numbering built in BuildPetscLocalMaps.
    ISLocalToGlobalMappingCreate(MPI_COMM_WORLD, this->ndof, cmesh.nodec, this->petsc_local_block_gids.data(), PETSC_COPY_VALUES,
                                 &this->petsc_block_lgmap_);
    ISLocalToGlobalMappingCreate(MPI_COMM_WORLD, 1, cmesh.nodec * this->ndof, this->petsc_local_var_gids.data(), PETSC_COPY_VALUES,
                                 &this->petsc_var_lgmap_);

    MatSetLocalToGlobalMapping(this->petsc_mat, this->petsc_block_lgmap_, this->petsc_block_lgmap_);
    VecSetLocalToGlobalMapping(this->petsc_b, this->petsc_var_lgmap_);
    VecSetLocalToGlobalMapping(this->petsc_x, this->petsc_var_lgmap_);

    VecZeroEntries(this->petsc_b);
    VecZeroEntries(this->petsc_x);

    // Block-oriented buffers for the optional blocked assembly path.
    this->petsc_block_cols_buf.resize(static_cast<size_t>(max_row_nnz));
    this->petsc_block_vals_buf.resize(static_cast<size_t>(max_row_nnz) * this->ndof * this->ndof);

    return;
}

void CrsMat::BuildKSPSolver() {
    KSPCreate(MPI_COMM_WORLD, &this->ksp);

    KSPSetOperators(this->ksp, this->petsc_mat, this->petsc_mat);

    KSPSetType(this->ksp, KSPFGMRES);

    PC pc;
    KSPGetPC(this->ksp, &pc);
    this->ConfigurePreconditioner(pc);

    KSPSetTolerances(this->ksp, 1.0e-10, 1.0e-15, 1.0e6, 1000);

    KSPSetFromOptions(this->ksp);

    return;
}

void CrsMat::InitPetscSolver() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->BuildLGMAP();

    this->BuildPetscMat();

    this->BuildKSPSolver();

    IS from_is = nullptr;
    IS to_is = nullptr;
    ISCreateGeneral(PETSC_COMM_SELF, cmesh.nodec * this->ndof, this->natural_var_gids.data(), PETSC_COPY_VALUES, &from_is);
    VecCreateSeq(PETSC_COMM_SELF, cmesh.nodec * this->ndof, &this->seq_x);
    ISCreateStride(PETSC_COMM_SELF, cmesh.nodec * this->ndof, 0, 1, &to_is);
    VecScatterCreate(this->petsc_x, from_is, this->seq_x, to_is, &this->scatter_to_all);
    ISDestroy(&from_is);
    ISDestroy(&to_is);

    return;
}
