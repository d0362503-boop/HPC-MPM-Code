#include <vector>
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "../../mesh.h"
#include "../../dataset.h"
#include "../../material_point.h"
#include "../../shape_function.h"
#include "../../cal_mat.h"
#include "../../mpi_data.h"
#include "../../solid/implicit/implicit_mpm_solid.h"
#include "fbar_projection.h"

using namespace std;
using namespace implicitmpm;

void cal_def_grad_bar (vector<array<array<double, 3>, 3>>& delta_def_grad,
                       const vector<double>& det_delta_def_grad) {

    int nenode;
    vector<int> ncm; 
    vector<double> sf;
    vector<array<double, 3>> dsf;

    vector<double> nJbar(nodec, 0.0e0);
    for (int m = 0; m < nelem; m++) {
        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = sp.coord[pid];

            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                nJbar[nid] += sfi * (sp.det_def_grad_bar[pid] * det_delta_def_grad[pid]) * sp.vol[pid];
            }
            pid = sp.idp2p[pid];
        }
    }

    NodeVarComm(nJbar, 0);

    for (int n = 0; n < nodec; n++) {
        if (sp.nvof[n] > mtol) {
            nJbar[n] /= sp.nvof[n];
        }
        else {
            nJbar[n] = 0.0e0;
        }
    }

    for (int m = 0; m < nelem; m++) {
        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = sp.coord[pid];

            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            double Jbar = 0.0e0;
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                Jbar += sfi * nJbar[nid];
            }

            double J = sp.det_def_grad_bar[pid] * det_delta_def_grad[pid];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    delta_def_grad[pid][i][j] *= cbrt((Jbar / J));
                }
            }

            std::array<std::array<double, 3>, 3> F_bar{};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        F_bar[i][j] += delta_def_grad[pid][i][k] * sp.def_grad_bar[pid][k][j];
                    }
                }
            }
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sp.def_grad_bar[pid][i][j] = F_bar[i][j];
                }
            }
            sp.det_def_grad_bar[pid] = Jbar;
            //det_mat3(F_bar, sp.det_def_grad_bar[pid]);

            pid = sp.idp2p[pid];
        }
    }

    return;
}
