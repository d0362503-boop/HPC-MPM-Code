#include <algorithm>
#include <mpi.h>
#include <utility>
#include <vector>

#include "module/cal_mat.h"
#include "module/DLB/mpm_dlb.h"
#include "module/dataset.h"
#include "module/interface_sdf.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::LumpedLagrangeMultiplier() {

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    double volp = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nlm_lump_local);
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
                        this->nlm_lump_local[nid] += sfi * volp * norm_grad_phi;
                    }
                }
            }
        }
    }

    VectorAssign(nodec, this->nlm_lump);
    this->nlm_lump = this->nlm_lump_local;

    NodeVarComm(this->nlm_lump, 0);

    return;
}

void MPMMPMMonolithicFSI::CorrectFluidPenetration() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    std::vector<double> final_phi(nodec, 0.0e0);
    for (int pid = 0; pid < this->solid_.num; pid++) {
        int cell;
        std::array<double, 3> xyp = this->solid_.coord[pid];
        if (!LocateLocalElement(xyp, cell)) continue;
        MakeSF(cell, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
        for (int ni = 0; ni < nenode; ni++) {
            const int nid = ncm[ni];
            const double sfi = sf[ni];
            this->solid_.StandardVarP2G(pid, nid, sfi, this->solid_.vol, final_phi);
        }
    }

    NodeVarComm(final_phi, 0);

    for (int n = 0; n < nodec; ++n) { final_phi[n] = std::clamp(final_phi[n] / nvol[n], 0.0, 1.0); }
    const interface_geometry::InterfaceSDF interface(std::move(final_phi));

    const auto &regions = mpm_dlb::CurrentRegions();
    std::vector<std::vector<int>> send_ids(nprocs);
    for (int pid = 0; pid < this->fluid_.num; ++pid) {
        auto &point = this->fluid_.coord[pid];
        int cell;
        if (LocateLocalElement(point, cell)) {
            std::array<double, 3> surface_point;
            if (interface.FindSurfacePoint(point, surface_point)) { point = surface_point; }
            continue;
        }
        std::array<int, 3> global_cell;
        if (!LocateGlobalElement(point, global_cell)) continue;
        for (int rank = 0; rank < nprocs; ++rank) {
            const auto &region = regions[rank];
            if (global_cell[0] >= region.elem_min[0] && global_cell[0] <= region.elem_max[0] &&
                global_cell[1] >= region.elem_min[1] && global_cell[1] <= region.elem_max[1] &&
                global_cell[2] >= region.elem_min[2] && global_cell[2] <= region.elem_max[2]) {
                send_ids[rank].push_back(pid);
                break;
            }
        }
    }

    ParticleCommunication query_comm;
    query_comm.nsp.resize(nprocs);
    query_comm.nrp.resize(nprocs);
    for (int rank = 0; rank < nprocs; ++rank) {
        query_comm.comm_ranks.push_back(rank);
        query_comm.nsp[rank] = send_ids[rank].size();
    }
    MPI_Alltoall(query_comm.nsp.data(), 1, MPI_INT, query_comm.nrp.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> query_particles;
    std::vector<std::array<double, 3>> send, receive;
    int receive_count = 0;
    for (int rank = 0; rank < nprocs; ++rank) {
        for (int pid : send_ids[rank]) {
            query_particles.push_back(pid);
            send.push_back(this->fluid_.coord[pid]);
        }
        receive_count += query_comm.nrp[rank];
    }
    receive.resize(receive_count);
    query_comm.PointVarSendrecv(send, receive, 1);
    for (auto &point : receive) {
        std::array<double, 3> surface_point;
        if (interface.FindSurfacePoint(point, surface_point)) { point = surface_point; }
    }

    // Reply in the same peer order, with the communication counts reversed.
    std::swap(query_comm.nsp, query_comm.nrp);
    query_comm.PointVarSendrecv(receive, send, 1);
    for (std::size_t i = 0; i < query_particles.size(); ++i) {
        this->fluid_.coord[query_particles[i]] = send[i];
    }

    return;
}
