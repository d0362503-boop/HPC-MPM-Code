#include <array>
#include <cmath>
#include <vector>

#include "../../cal_mat.h"
#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "explicit_mpm_solid.h"

using namespace explicitmpm;

void ExplicitSolidMPM::ComputeDefGradBar(const std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                                         const std::vector<double> &det_delta_def_grad) {
    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    // --- Step 1: project uncorrected new Jacobian (weighted by particle volume) to control points ---
    std::vector<double> nJbar(nodec, 0.0e0);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            const double J_new = this->det_def_grad_bar[pid] * det_delta_def_grad[pid];
            const double J_new_vol = J_new * this->vol[pid];

            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                const int nid = ncm[ni];
                const double sfi = sf[ni];
                nJbar[nid] += sfi * J_new_vol;
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(nJbar, 0);

    this->CutOffSmallNodalVar(nJbar, this->nvof, {0});

    // --- Step 2: interpolate corrected Jacobian and update F-bar deformation gradient ---
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            double Jbar = 0.0e0;
            for (int ni = 0; ni < nenode; ni++) {
                const int nid = ncm[ni];
                const double sfi = sf[ni];
                Jbar += sfi * nJbar[nid];
            }

            const double J = this->det_def_grad_bar[pid] * det_delta_def_grad[pid];
            const double scale = std::cbrt(Jbar / J);
            auto dF = delta_def_grad[pid];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) { dF[i][j] *= scale; }
            }

            std::array<std::array<double, 3>, 3> F_bar{};
            const auto F_bar_old = this->def_grad_bar[pid];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) { F_bar[i][j] += dF[i][k] * F_bar_old[k][j]; }
                }
            }
            this->def_grad_bar[pid] = F_bar;
            this->det_def_grad_bar[pid] = Jbar;

            pid = this->idp2p[pid];
        }
    }

    return;
}
