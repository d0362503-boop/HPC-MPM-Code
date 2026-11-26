#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"
#include "module/solver/crsmat.h"

using namespace stabilizedfem;

void StabilizedFEM::ConfigurePreconditioner(CrsMat &mat, PC pc) {

    if (mat.ndof == 1) {
        this->MaterialPoint::ConfigurePreconditioner(mat, pc);
        return;
    }

    const PetscInt velocity_fields[] = {0, 1, 2};
    const PetscInt pressure_field = 3;

    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetBlockSize(pc, mat.ndof);
    PCFieldSplitSetFields(pc, "velocity", 3, velocity_fields, velocity_fields);
    PCFieldSplitSetFields(pc, "pressure", 1, &pressure_field, &pressure_field);
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);

    PetscOptionsSetValue(nullptr, "-fieldsplit_velocity_ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-fieldsplit_velocity_pc_type", "hypre");
    PetscOptionsSetValue(nullptr, "-fieldsplit_velocity_pc_hypre_type", "boomeramg");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_type", "hypre");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_type", "boomeramg");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_boomeramg_coarsen_type", "hmis");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_boomeramg_interp_type", "ext+i");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_boomeramg_max_iter", "2");

    return;
}

void StabilizedFEM::BeginNewtonStep() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    VectorAssign(cmesh.nodec * 3, this->nvel_k);
    VectorAssign(cmesh.nodec, this->npres_k);
    this->r0r = 0.0e0;

    std::copy(this->nvel.begin(), this->nvel.end(), this->nvel_k.begin());
    std::copy(this->npres.begin(), this->npres.end(), this->npres_k.begin());

    return;
}

std::vector<double> StabilizedFEM::ComputeAdvectionVel() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    std::vector<double> adv_vel(cmesh.nodec * 3);
    for (int n = 0; n < cmesh.nodec * 3; n++) {
        // --- v^{n+α_f} from level n and the n+1 iterate ---
        adv_vel[n] = (1.0e0 - this->integrator_.alpha_f) * this->nvel[n] //
                     + this->integrator_.alpha_f * this->nvel_k[n];
    }

    return adv_vel;
}

void StabilizedFEM::AssembleNewtonSystem(const TimeState &ts) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->BCNRSet(ts);

    std::vector<double> adv_vel = this->ComputeAdvectionVel();

    this->naccel_k = this->ComputeNodeAccelFromVel();

    this->MakeNSStabCoeff(adv_vel); // ---- Stabilized coefficient ----

    VectorAssign(this->NS_.nmata, this->NS_.amat);
    VectorAssign(cmesh.nodec * 4, this->NS_.b_rhs);
    VectorAssign(cmesh.nodec * 4, this->NS_.x_lhs);
    this->AssembleSystem(this->NS_, adv_vel, ts);

    return;
}

void StabilizedFEM::SolveNS(const TimeState &ts) {

    SolveNewton(*this, ts, this->max_NR_it);

    return;
}

bool StabilizedFEM::CheckNewtonConvergence(int NR_it, int ksp_its) {

    return this->NS_.CheckNRConvergence(NR_it, this->max_NR_it, ksp_its, this->r0r);
}

int StabilizedFEM::SolveNewtonCorrection(int NR_it, const TimeState &ts) { return this->NS_.SolveSystem(NR_it, ts.istep); }

void StabilizedFEM::EndNewtonStep() {

    // --- Memory clear ---
    std::vector<double>().swap(this->NS_.adiag);
    std::vector<double>().swap(this->NS_.amat);
    std::vector<double>().swap(this->NS_.b_rhs);
    std::vector<double>().swap(this->NS_.x_lhs);

    return;
}

void StabilizedFEM::MakeNSStabCoeff(const std::vector<double> &adv_vel) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    this->tau_supg.resize(mesh.nelem);
    this->tau_lsic.resize(mesh.nelem);
    for (int m = 0; m < mesh.nelem; m++) {
        double rnue = this->rmue[m] / this->rhoe[m];

        const std::array<double, 3> xye = ElementCenter(mesh, m);

        MakeSF(this->grid_, m, xye, ncm, nenode, sf, dsf);

        double uu = 0.0e0, vv = 0.0e0, ww = 0.0e0;
        for (int ni = 0; ni < nenode; ni++) {
            int nid = ncm[ni];
            uu += adv_vel[nid + cmesh.nuc] / double(nenode);
            vv += adv_vel[nid + cmesh.nvc] / double(nenode);
            ww += adv_vel[nid + cmesh.nwc] / double(nenode);
        }
        double uvw = uu * uu + vv * vv + ww * ww;

        double udsf = 0.0e0;
        for (int ni = 0; ni < nenode; ni++) {
            double dsfi1 = dsf[ni][0];
            double dsfi2 = dsf[ni][1];
            double dsfi3 = dsf[ni][2];
            udsf += std::abs(uu * dsfi1 + vv * dsfi2 + ww * dsfi3);
        }

        double he;
        if (udsf < params.mtol) {
            he = mesh.dxy[0] * std::cbrt((6.0e0 / M_PI));
        } else {
            he = 2.0e0 * std::sqrt(uvw) / udsf;
        }

        double t1 = (2.0e0 * params.dti) * (2.0e0 * params.dti);
        double t2 = 4.0e0 * uvw / (he * he);
        double t3 = std::pow((4.0e0 * rnue / (he * he)), 2);
        // --- τ_SUPG = τ_PSPG = [(2/Δt)^2 + (2|v|/h)^2 + (4ν/h^2)^2]^{-1/2} ---
        this->tau_supg[m] = 1.0e0 / std::sqrt(t1 + t2 + t3);

        // --- τ_LSIC = (h/2)|v| z(Re_h), z = Re_h/3 for Re_h <= 3, else 1 ---
        if (he < params.mtol) {
            this->tau_lsic[m] = 0.0e0;
        } else {
            double reu = std::sqrt(uvw) * he / (rnue * 2.0e0);
            if (reu > 3.0e0) {
                this->tau_lsic[m] = std::sqrt(uvw) * he / 2.0e0;
            } else {
                this->tau_lsic[m] = std::sqrt(uvw) * he * reu / 6.0e0;
            }
        }
    }

    return;
}

void StabilizedFEM::AssembleSystem(CrsMat &mat, const std::vector<double> &adv_vel, const TimeState &ts) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    double fx = params.bb[0] * ts.facl;
    double fy = params.bb[1] * ts.facl;
    double fz = params.bb[2] * ts.facl;

    const std::vector<int> offsets = this->RHSOffsets();

    const double am = this->integrator_.alpha_m;
    const double af = this->integrator_.alpha_f;
    const double af0 = 1.0e0 - af;
    const double am0 = 1.0e0 - am;
    const double gamma = this->integrator_.gamma_nb;
    const double c_m = am / (gamma * params.dt); // d a^{n+α_m} / d v^{n+1}

    // --- Convected state at its time levels; adv_vel is the convecting velocity
    //     (equal to nvel_af on a fixed Eulerian grid, v - v_mesh under ALE) ---
    std::vector<double> nvel_af(cmesh.nodec * 3), naccel_am(cmesh.nodec * 3);
    for (int n = 0; n < cmesh.nodec * 3; n++) {
        naccel_am[n] = am0 * this->naccel[n] + am * this->naccel_k[n];
        nvel_af[n] = af0 * this->nvel[n] + af * this->nvel_k[n];
    }

    std::vector<double> npres_af(cmesh.nodec);
    for (int n = 0; n < cmesh.nodec; n++) { //
        npres_af[n] = af0 * this->npres[n] + af * this->npres_k[n];
    }

    double weight = (mesh.dxy[0] * mesh.dxy[1] * mesh.dxy[2]) / (params.npxye[0] * params.npxye[1] * params.npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    for (int m = 0; m < mesh.nelem; m++) {
        double t1 = this->tau_supg[m];
        double t2 = this->tau_lsic[m];

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(this->grid_, m, xyp, ncm, nenode, sf, dsf);

            double phi_k = 0.0e0;
            std::array<double, 3> adv_vel_k{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                // --- Liquid volume fraction ---
                phi_k += sfi * this->nphi[nid];
                // --- Advection velocity ---
                adv_vel_k[0] += sfi * adv_vel[nid + cmesh.nuc];
                adv_vel_k[1] += sfi * adv_vel[nid + cmesh.nvc];
                adv_vel_k[2] += sfi * adv_vel[nid + cmesh.nwc];
            }

            double rhof = phi_k * this->rhol + (1.0e0 - phi_k) * this->rhog;
            double rmuf = phi_k * this->rmul + (1.0e0 - phi_k) * this->rmug;

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                double ust_ni = (dsfi1 * adv_vel_k[0] + dsfi2 * adv_vel_k[1] + dsfi3 * adv_vel_k[2]); // v·∇N_i
                // --- Force (Galerkin+SUPG): ∫ (N_i + τ v·∇N_i) ρ f ---
                double fg = (sfi + ust_ni * t1) * weight * rhof;
                // --- Force (PSPG): ∫ (τ/ρ) ∇N_i · ρ f ---
                double fp = (dsfi1 * fx + dsfi2 * fy + dsfi3 * fz) * weight * t1;

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];
                    double ust_nj = (dsfj1 * adv_vel_k[0] + dsfj2 * adv_vel_k[1] + dsfj3 * adv_vel_k[2]); // v·∇N_j

                    int ida = mat.FindIndex(nid, njd, ncol);

                    // Entries below are unscaled operators (i: test node, j: trial node, k: direction).
                    // Momentum rows use the SUPG test function N_i + τ v·∇N_i; continuity rows
                    // use N_i (Galerkin) plus (τ/ρ)∇N_i (PSPG).

                    // --- Mass (Galerkin+SUPG): ∫ (N_i + τ v·∇N_i) ρ N_j, acts on a^{n+α_m} ---
                    double emd = weight * (sfi + ust_ni * t1) * sfj * rhof;
                    // --- Mass (PSPG): ∫ τ ∂N_i/∂x_k N_j, acts on a_k^{n+α_m} ---
                    double egtu = weight * dsfi1 * sfj * t1;
                    double egtv = weight * dsfi2 * sfj * t1;
                    double egtw = weight * dsfi3 * sfj * t1;
                    // --- Advection (Galerkin+SUPG): ∫ (N_i + τ v·∇N_i) ρ (v·∇N_j) ---
                    double akg = weight * (sfi + ust_ni * t1) * ust_nj * rhof;
                    // --- Advection (PSPG): ∫ τ ∂N_i/∂x_k (v·∇N_j), acts on v_k ---
                    double akpu = weight * dsfi1 * ust_nj * t1;
                    double akpv = weight * dsfi2 * ust_nj * t1;
                    double akpw = weight * dsfi3 * ust_nj * t1;
                    // --- Viscous (Galerkin): ∫ ε(N_i e_k) : 2μ ε(N_j e_l) ---
                    //     kk: μ(∇N_i·∇N_j + ∂N_i/∂x_k ∂N_j/∂x_k), kl: μ ∂N_i/∂x_l ∂N_j/∂x_k
                    double su = weight * (2.0e0 * dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * rmuf;
                    double sv = weight * (dsfi1 * dsfj1 + 2.0e0 * dsfi2 * dsfj2 + dsfi3 * dsfj3) * rmuf;
                    double sw = weight * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + 2.0e0 * dsfi3 * dsfj3) * rmuf;
                    double suv = weight * (dsfi2 * dsfj1) * rmuf;
                    double suw = weight * (dsfi3 * dsfj1) * rmuf;
                    double svu = weight * (dsfi1 * dsfj2) * rmuf;
                    double svw = weight * (dsfi3 * dsfj2) * rmuf;
                    double swu = weight * (dsfi1 * dsfj3) * rmuf;
                    double swv = weight * (dsfi2 * dsfj3) * rmuf;
                    // --- Pressure gradient: ∫ [τ (v·∇N_i) ∂N_j/∂x_k − ∂N_i/∂x_k N_j] (SUPG − Galerkin (∇·w)p) ---
                    double esgu = weight * (ust_ni * dsfj1 * t1 - dsfi1 * sfj);
                    double esgv = weight * (ust_ni * dsfj2 * t1 - dsfi2 * sfj);
                    double esgw = weight * (ust_ni * dsfj3 * t1 - dsfi3 * sfj);
                    // --- Pressure (PSPG): ∫ (τ/ρ) ∇N_i·∇N_j ---
                    double elt = weight * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) / rhof * t1;
                    // --- Continuity (Galerkin): ∫ N_i ∂N_j/∂x_k, i.e. q ∇·v ---
                    double Cowu = weight * sfi * dsfj1;
                    double Cowv = weight * sfi * dsfj2;
                    double Coww = weight * sfi * dsfj3;
                    // --- LSIC (grad-div): ∫ τ_LSIC ρ ∂N_i/∂x_k ∂N_j/∂x_l, i.e. (∇·w)(∇·v) ---
                    double scu = weight * dsfi1 * dsfj1 * rhof * t2;
                    double scv = weight * dsfi2 * dsfj2 * rhof * t2;
                    double scw = weight * dsfi3 * dsfj3 * rhof * t2;
                    double scuv = weight * dsfi1 * dsfj2 * rhof * t2;
                    double scuw = weight * dsfi1 * dsfj3 * rhof * t2;
                    double scvu = weight * dsfi2 * dsfj1 * rhof * t2;
                    double scvw = weight * dsfi2 * dsfj3 * rhof * t2;
                    double scwu = weight * dsfi3 * dsfj1 * rhof * t2;
                    double scwv = weight * dsfi3 * dsfj2 * rhof * t2;

                    // --- Picard tangent dR/d(v^{n+1}, p^{n+1}); dτ/dv and (δv·∇)v are omitted ---
                    //     acceleration terms × c_m, all velocity and pressure terms × α_f.
                    //     blocks[4*r + c]: row r, column c in (u, v, w, p).
                    mat.amat[ida + mat.block_id[this->blocks[0]]] += c_m * emd + af * (akg + su + scu); // u-u
                    mat.amat[ida + mat.block_id[this->blocks[1]]] += af * (suv + scuv);                // u-v
                    mat.amat[ida + mat.block_id[this->blocks[2]]] += af * (suw + scuw);                // u-w
                    mat.amat[ida + mat.block_id[this->blocks[3]]] += af * esgu;                        // u-p
                    mat.amat[ida + mat.block_id[this->blocks[4]]] += af * (svu + scvu);                // v-u
                    mat.amat[ida + mat.block_id[this->blocks[5]]] += c_m * emd + af * (akg + sv + scv); // v-v
                    mat.amat[ida + mat.block_id[this->blocks[6]]] += af * (svw + scvw);                // v-w
                    mat.amat[ida + mat.block_id[this->blocks[7]]] += af * esgv;                        // v-p
                    mat.amat[ida + mat.block_id[this->blocks[8]]] += af * (swu + scwu);                // w-u
                    mat.amat[ida + mat.block_id[this->blocks[9]]] += af * (swv + scwv);                // w-v
                    mat.amat[ida + mat.block_id[this->blocks[10]]] += c_m * emd + af * (akg + sw + scw); // w-w
                    mat.amat[ida + mat.block_id[this->blocks[11]]] += af * esgw;                        // w-p
                    mat.amat[ida + mat.block_id[this->blocks[12]]] += c_m * egtu + af * (Cowu + akpu);  // p-u
                    mat.amat[ida + mat.block_id[this->blocks[13]]] += c_m * egtv + af * (Cowv + akpv);  // p-v
                    mat.amat[ida + mat.block_id[this->blocks[14]]] += c_m * egtw + af * (Coww + akpw);  // p-w
                    mat.amat[ida + mat.block_id[this->blocks[15]]] += af * elt;                         // p-p

                    // --- Residual RHS = -R: the same unscaled entries applied to
                    //     a^{n+α_m}, v^{n+α_f}, p^{n+α_f}, plus the body force below ---

                    const double u_af = nvel_af[njd + cmesh.nuc];
                    const double v_af = nvel_af[njd + cmesh.nvc];
                    const double w_af = nvel_af[njd + cmesh.nwc];
                    const double p_af = npres_af[njd];

                    const double au_am = naccel_am[njd + cmesh.nuc];
                    const double av_am = naccel_am[njd + cmesh.nvc];
                    const double aw_am = naccel_am[njd + cmesh.nwc];

                    // x-momentum
                    mat.b_rhs[nid + offsets[0]] += -(akg + su + scu) * u_af //
                                                   - (suv + scuv) * v_af    //
                                                   - (suw + scuw) * w_af    //
                                                   - esgu * p_af - emd * au_am;

                    // y-momentum
                    mat.b_rhs[nid + offsets[1]] += -(akg + sv + scv) * v_af //
                                                   - (svu + scvu) * u_af    //
                                                   - (svw + scvw) * w_af    //
                                                   - esgv * p_af - emd * av_am;

                    // z-momentum
                    mat.b_rhs[nid + offsets[2]] += -(akg + sw + scw) * w_af //
                                                   - (swu + scwu) * u_af    //
                                                   - (swv + scwv) * v_af    //
                                                   - esgw * p_af - emd * aw_am;

                    // continuity / pressure
                    mat.b_rhs[nid + offsets[3]] += -(Cowu + akpu) * u_af  //
                                                   - (Cowv + akpv) * v_af //
                                                   - (Coww + akpw) * w_af //
                                                   - elt * p_af - (egtu * au_am + egtv * av_am + egtw * aw_am);
                }
                mat.b_rhs[nid + offsets[0]] += fg * fx;
                mat.b_rhs[nid + offsets[1]] += fg * fy;
                mat.b_rhs[nid + offsets[2]] += fg * fz;
                mat.b_rhs[nid + offsets[3]] += fp;
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->grid_.part, mat.b_rhs, offsets);

    return;
}

void StabilizedFEM::UpdateNRIncrement() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    for (int n = 0; n < cmesh.nodec * 3; n++) { this->nvel_k[n] += this->NS_.x_lhs[n]; }
    for (int n = 0; n < cmesh.nodec; n++) { this->npres_k[n] += this->NS_.x_lhs[n + cmesh.npc]; }

    return;
}

void StabilizedFEM::CommitNodalSolution() {

    this->naccel = this->ComputeNodeAccelFromVel();
    std::copy(this->nvel_k.begin(), this->nvel_k.end(), this->nvel.begin());
    std::copy(this->npres_k.begin(), this->npres_k.end(), this->npres.begin());

    return;
}
