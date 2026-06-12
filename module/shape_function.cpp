#include <vector>
#include <array>
#include "mesh.h"
#include "dataset.h"
#include "material_point.h"
#include "cal_mat.h"

void MakSf(const int m, const std::array<double, 3>& xyg, const std::vector<int>& idimc, const std::vector<int>& xynodec,
           std::vector<int>& nc, int& nenode, std::vector<double>& sf, std::vector<std::array<double, 3>>& dsf) {

    const int max_nenode = (idimc[0]+1) * (idimc[1]+1) * (idimc[2]+1);
    int nxc = xynodec[0];
    int nyc = xynodec[1];
    int nzc = xynodec[2];
    int nex = xyelem[0];
    int ney = xyelem[1];
    int nez = xyelem[2];

    int ize = m / (nex * ney);
    int iye = (m - ize * (nex * ney)) / nex;
    int ixe = m - ize * (nex * ney) - iye * nex;

    nc.resize(max_nenode); sf.resize(max_nenode);
    dsf.resize(max_nenode);

    int ixyn[3];
    int id = 0;
    for (int k = 0; k <= idimc[2]; k++) {
        ixyn[2] = ize + k;
        for (int j = 0; j <= idimc[1]; j++) {
            ixyn[1] = iye + j;
            for (int i = 0; i <= idimc[0]; i++) {
                ixyn[0] = ixe + i;
                nc[id] = ixyn[0] + nxc * ixyn[1] + nxc * nyc * ixyn[2];

                std::array<double, 3> sfxy, dsfxy;
                for (int ii = 0; ii < 3; ++ii) {
                    double dxi = dxy[ii];
                    // ----- Linear shape function -----
                    if (idimc[ii] == 1) {
                        double xy = xymin[ii] + dxi * double(ixyn[ii]);
                        double dxnpi = xyg[ii] - xy;
                        double xi = abs(dxnpi) / dxi;
                        sfxy[ii] = 1.0e0 - xi;
                        if (dxnpi < 0.0e0) {
                            dsfxy[ii] = 1.0e0 / dxi;
                        }
                        else {
                            dsfxy[ii] = -1.0e0 / dxi;
                        }
                    }
                    // ----- Second order B-spline shape function -----
                    else if (idimc[ii] == 2) {
                        int ncxy = aelemmin[ii] + ixyn[ii];
                        double xy = xymin[ii] + dxi * double(ixyn[ii] - 2);
                        double dxnpi = xyg[ii] - xy;
                        double dxh = dxnpi / dxi;
                        double s1, s2, s3, d1, d2;
                        if (ncxy == 0) {
                            dxh -= 2.0e0;
                            s1 = 1.0e0; s2 = -2.0e0; s3 = 1.0e0;
                            d1 = 2.0e0; d2 = -2.0e0;
                        }
                        else if (ncxy == 1) {
                            dxh -= 1.0e0;
                            if (dxh <= 1.0e0) {
                                s1 = -1.5e0; s2 = 2.0e0; s3 = 0.0e0;
                                d1 = -3.0e0; d2 = 2.0e0;
                            }
                            else {
                                s1 = 0.5e0; s2 = -2.0e0; s3 = 2.0e0;
                                d1 = 1.0e0; d2 = -2.0e0;
                            }
                        }
                        else if (ncxy == xynodecw[ii] - 1) {
                            s1 = 1.0e0; s2 = 0.0e0; s3 = 0.0e0;
                            d1 = 2.0e0; d2 = 0.0e0;
                        }
                        else if (ncxy == xynodecw[ii] - 2) {
                            if (dxh <= 1.0e0) {
                                s1 = 0.5e0; s2 = 0.0e0; s3 = 0.0e0;
                                d1 = 1.0e0; d2 = 0.0e0;
                            }
                            else {
                                s1 = -1.5e0; s2 = 4.0e0; s3 = -2.0e0;
                                d1 = -3.0e0; d2 = 4.0e0;
                            }
                        }
                        else {
                            if (dxh <= 1.0e0) {
                                s1 = 0.5e0; s2 = 0.0e0; s3 = 0.0e0;
                                d1 = 1.0e0; d2 = 0.0e0;
                            }
                            else if (dxh <= 2.0e0) {
                                s1 = -1.0e0; s2 = 3.0e0; s3 = -1.5e0;
                                d1 = -2.0e0; d2 = 3.0e0;
                            }
                            else {
                                s1 = 0.5e0; s2 = -3.0e0; s3 = 4.5e0;
                                d1 = 1.0e0; d2 = -3.0e0;
                            }
                        }
                        sfxy[ii] = (s1 * dxh + s2) * dxh + s3;
                        dsfxy[ii] = (d1 * dxh + d2) / dxi;
                    }
                }
                sf[id] = sfxy[0] * sfxy[1] * sfxy[2];
                dsf[id][0] = dsfxy[0] * sfxy[1] * sfxy[2];
                dsf[id][1] = sfxy[0] * dsfxy[1] * sfxy[2];
                dsf[id][2] = sfxy[0] * sfxy[1] * dsfxy[2];

                id++;
            }
        }
    }
    nenode = id;

    return;
}
