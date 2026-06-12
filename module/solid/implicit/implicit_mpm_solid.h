#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../solver/crsmat.h"
#include "../solid_material_point.h"

namespace implicitmpm {

class SolidMaterialPoint : public SolidMaterialPointBase {
  public:
    CrsMat SM_;

    SolidMaterialPoint() {
        this->ode_order = 2;
        this->cm_.implicit_flag = true;
        this->SM_.ndof = 3;
        this->SM_.use_petsc = true;
        this->SM_.FEM_flag = false;
        this->SM_.amg_rebuild_freq = 1; // keep solid AMG alive across this 10-step window unless iterations deteriorate
        this->SM_.owner_ = this;
    }

    void SetTracForce() {
        const double t0 = 1.0e0 * dxy[1] * dxy[2] / (npxye[1] * npxye[2]);
        double t = 0.0e0;
        if (real_time < 5.0e-3) { t = -t0; }

        for (int ip = 0; ip < this->num; ip++) {
            if (this->surf_point[ip] == 1) { this->trac_force[ip][0] = t; }
        }

        return;
    };

    void Particle2Node() override;

    void Node2Particle() override;

    void SolveSolid() override;

    void UpdateDefGrad(int pid, int nenode, double af_coeff, const std::vector<int> &ncm, const std::vector<double> &sf,
                       const std::vector<std::array<double, 3>> &dsf,
                       std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                       std::vector<std::array<std::array<double, 3>, 3>> &def_grad);

    std::vector<std::array<double, 6>> InitializeNRStress();

    void BCNRSet() override {
        this->ubc.BCSetDt(nuc, this->ndispl);
        this->vbc.BCSetDt(nvc, this->ndispl);
        this->wbc.BCSetDt(nwc, this->ndispl);
        this->rigid_bc.BCSetDt(nuc, this->ndispl);
        this->rigid_bc.BCSetDt(nvc, this->ndispl);
        this->rigid_bc.BCSetDt(nwc, this->ndispl);

        return;
    };

    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(nuc, rr);
        this->vbc.BCSetZero(nvc, rr);
        this->wbc.BCSetZero(nwc, rr);
        this->rigid_bc.BCSetZero(nuc, rr);
        this->rigid_bc.BCSetZero(nvc, rr);
        this->rigid_bc.BCSetZero(nwc, rr);

        return;
    };

    void BuildPetscBCList(CrsMat &mat) override {
        mat.AddBCComponent(this->ubc, nuc);
        mat.AddBCComponent(this->vbc, nvc);
        mat.AddBCComponent(this->wbc, nwc);
        mat.AddBCComponent(this->rigid_bc, nuc);
        mat.AddBCComponent(this->rigid_bc, nvc);
        mat.AddBCComponent(this->rigid_bc, nwc);

        return;
    };

    void AssembleSystem(const std::vector<double> &naccel_k, const std::vector<double> &nvel_k,
                        std::vector<std::array<double, 6>> &stress_k);

    auto ComputeTangentModulus(int pid, int ni, int nj, const std::vector<std::array<double, 3>> &dsf,
                               const std::array<double, 6> &sts_af);

    void UpdateNRIncrement() override;
};

// --- Solid Point object ---
inline SolidMaterialPoint sp;
} // namespace implicitmpm
