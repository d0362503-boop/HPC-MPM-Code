#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <utility>
#include <vector>

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/interface_sdf.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace mpm_mpm_monolithic_fsi;

template <size_t vertex_count>
std::vector<std::array<std::array<double, 3>, vertex_count>>
MPMMPMMonolithicFSI::GatherInterfaceGeometry(const std::vector<std::array<std::array<double, 3>, vertex_count>> &local_geometry)
    const {
    std::vector<double> packed;
    for (const auto &element : local_geometry) {
        for (const auto &vertex : element) { packed.insert(packed.end(), vertex.begin(), vertex.end()); }
    }

    const int local_count = packed.size();
    std::vector<int> counts(nprocs), offsets(nprocs);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int rank = 1; rank < nprocs; ++rank) { offsets[rank] = offsets[rank - 1] + counts[rank - 1]; }

    std::vector<double> gathered(offsets.back() + counts.back());
    MPI_Allgatherv(packed.data(), local_count, MPI_DOUBLE, gathered.data(), counts.data(), offsets.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);

    std::vector<std::array<std::array<double, 3>, vertex_count>> geometry(gathered.size() / (3 * vertex_count));
    size_t index = 0;
    for (auto &element : geometry) {
        for (auto &vertex : element) {
            for (double &coordinate : vertex) { coordinate = gathered[index++]; }
        }
    }

    return geometry;
}

std::array<double, 3> MPMMPMMonolithicFSI::ContourIntersection(const std::array<std::array<double, 3>, 4> &vertices,
                                                               const std::array<double, 4> &level, int i, int j) const {

    const double fraction = level[i] / (level[i] - level[j]);
    std::array<double, 3> point;
    for (int d = 0; d < 3; d++) { point[d] = vertices[i][d] + fraction * (vertices[j][d] - vertices[i][d]); }

    return point;
}

void MPMMPMMonolithicFSI::AppendInterfaceTriangle(std::array<std::array<double, 3>, 3> triangle,
                                                  const std::array<double, 3> &outward,
                                                  std::vector<std::array<std::array<double, 3>, 3>> &surface) const {

    const std::array<double, 3> normal =
        CrossVec3(DifferenceVec3(triangle[1], triangle[0]), DifferenceVec3(triangle[2], triangle[0]));

    // A contour through a vertex can produce a zero-area triangle.
    if (DotVec3(normal, normal) == 0.0) return;
    if (DotVec3(normal, outward) < 0.0) { std::swap(triangle[1], triangle[2]); }

    surface.push_back(triangle);

    return;
}

void MPMMPMMonolithicFSI::AppendSolidContour(const std::array<std::array<double, 3>, 4> &vertices,
                                             const std::array<double, 4> &level,
                                             std::vector<std::array<std::array<double, 3>, 3>> &surface) const {

    std::vector<int> inside, outside;
    for (int i = 0; i < 4; i++) {
        if (level[i] < 0.0) {
            inside.push_back(i);
        } else {
            outside.push_back(i);
        }
    }
    if (inside.empty() || outside.empty()) return;

    std::array<double, 3> outward{};
    for (int d = 0; d < 3; d++) {
        for (int i : outside) { outward[d] += vertices[i][d] / outside.size(); }
        for (int i : inside) { outward[d] -= vertices[i][d] / inside.size(); }
    }

    if (inside.size() == 3) { inside.swap(outside); }

    if (inside.size() == 1) {
        std::array<std::array<double, 3>, 3> triangle;
        for (int i = 0; i < 3; i++) { triangle[i] = this->ContourIntersection(vertices, level, inside[0], outside[i]); }

        this->AppendInterfaceTriangle(triangle, outward, surface);
    } else {
        const std::array<double, 3> ac = this->ContourIntersection(vertices, level, inside[0], outside[0]);
        const std::array<double, 3> ad = this->ContourIntersection(vertices, level, inside[0], outside[1]);
        const std::array<double, 3> bc = this->ContourIntersection(vertices, level, inside[1], outside[0]);
        const std::array<double, 3> bd = this->ContourIntersection(vertices, level, inside[1], outside[1]);

        this->AppendInterfaceTriangle({ac, ad, bc}, outward, surface);
        this->AppendInterfaceTriangle({ad, bd, bc}, outward, surface);
    }

    return;
}

std::vector<std::array<std::array<double, 3>, 3>>
MPMMPMMonolithicFSI::BuildSolidInterface(const std::vector<double> &solid_phi) const {

    // Tetrahedra follow the nc corner order.
    const std::array<std::array<int, 4>, 6> tetrahedra{
        {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}}};
    std::vector<std::array<std::array<double, 3>, 3>> surface;
    std::vector<double> node_phi(node, 0.0e0);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    for (int m = 0; m < nelem; m++) {
        for (int n = 0; n < 8; n++) {
            const int id = nc[m][n];
            MakeSF(m, xyn[id], idimc, xynodec, ncm, nenode, sf, dsf);

            // Interpolate relative to the contour.
            node_phi[id] = 0.0e0;
            for (int ni = 0; ni < nenode; ni++) {
                const int nid = ncm[ni];
                const double sfi = sf[ni];
                node_phi[id] += sfi * (solid_phi[nid] - 0.5);
            }
            node_phi[id] += 0.5;
        }
    }

    for (int m = 0; m < nelem; m++) {
        double minimum_phi = 1.0, maximum_phi = 0.0;
        for (int id : nc[m]) {
            minimum_phi = std::min(minimum_phi, node_phi[id]);
            maximum_phi = std::max(maximum_phi, node_phi[id]);
        }
        if (minimum_phi > 0.5 || maximum_phi <= 0.5) continue;

        for (const std::array<int, 4> &tetrahedron : tetrahedra) {
            std::array<std::array<double, 3>, 4> positions;
            std::array<double, 4> level;
            for (int corner = 0; corner < 4; corner++) {
                const int id = nc[m][tetrahedron[corner]];
                positions[corner] = xyn[id];
                level[corner] = 0.5 - node_phi[id];
            }

            this->AppendSolidContour(positions, level, surface);
        }
    }

    return this->GatherInterfaceGeometry(surface);
}

std::vector<std::array<std::array<double, 3>, 2>> MPMMPMMonolithicFSI::BuildFluidDomains() const {

    std::vector<std::array<std::array<double, 3>, 2>> domains;
    const double reference_volume = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    for (int cell = 0; cell < nelem; ++cell) {
        double solid_support = 0.0;
        for (int node_id : ncc[cell]) { solid_support = std::max(solid_support, this->solid_.nphi[node_id]); }
        if (solid_support == 0.0) continue;

        int pid = this->fluid_.idepf[cell];
        while (pid != -1) {
            const double volume_scale = std::cbrt(this->fluid_.vol[pid] / reference_volume);
            std::array<std::array<double, 3>, 2> domain;

            // ponytail: preserve volume and grid aspect ratio; fluid domain deformation is not tracked.
            for (int d = 0; d < 3; ++d) {
                const double half_width = 0.5 * dxy[d] / npxye[d] * volume_scale;
                domain[0][d] = this->fluid_.coord[pid][d] - half_width;
                domain[1][d] = this->fluid_.coord[pid][d] + half_width;
            }

            domains.push_back(domain);

            pid = this->fluid_.idp2p[pid];
        }
    }

    return this->GatherInterfaceGeometry(domains);
}

std::map<std::array<int, 3>, std::vector<std::size_t>>
MPMMPMMonolithicFSI::IndexFluidDomains(const std::vector<std::array<std::array<double, 3>, 2>> &domains) const {

    std::map<std::array<int, 3>, std::vector<std::size_t>> bins;
    for (std::size_t domain_id = 0; domain_id < domains.size(); domain_id++) {
        std::array<int, 3> lower, upper;
        for (int d = 0; d < 3; d++) {
            const double margin = 1.5 * dxy[d];
            lower[d] = int(std::floor((domains[domain_id][0][d] - margin) / dxy[d]));
            upper[d] = int(std::floor((domains[domain_id][1][d] + margin) / dxy[d]));
        }

        for (int k = lower[2]; k <= upper[2]; k++) {
            for (int j = lower[1]; j <= upper[1]; j++) {
                for (int i = lower[0]; i <= upper[0]; i++) { bins[{i, j, k}].push_back(domain_id); }
            }
        }
    }

    return bins;
}

double MPMMPMMonolithicFSI::WetInterfaceAreaDensity(const std::array<double, 3> &point,
                                                    const interface_geometry::InterfaceSDF &interface,
                                                    const std::vector<std::array<std::array<double, 3>, 2>> &domains,
                                                    const std::map<std::array<int, 3>, std::vector<std::size_t>> &bins) const {

    std::array<double, 3> surface_point, normal;
    double signed_distance;
    if (!interface.FindSurfacePoint(point, surface_point, signed_distance, normal)) return 0.0;

    std::array<double, 3> normal_spacing;
    for (int d = 0; d < 3; d++) { normal_spacing[d] = normal[d] * dxy[d]; }

    const double band_width = NormVec3(normal_spacing);
    if (std::abs(signed_distance) >= band_width) return 0.0;

    const std::array<int, 3> cell{int(std::floor(surface_point[0] / dxy[0])), int(std::floor(surface_point[1] / dxy[1])),
                                  int(std::floor(surface_point[2] / dxy[2]))};
    const auto fluid_bin = bins.find(cell);
    if (fluid_bin == bins.end()) return 0.0;

    for (std::size_t domain_id : fluid_bin->second) {
        std::array<double, 3> separation{};
        for (int d = 0; d < 3; d++) {
            const double nearest = std::clamp(surface_point[d], domains[domain_id][0][d], domains[domain_id][1][d]);
            separation[d] = (surface_point[d] - nearest) / (1.5 * dxy[d]);
        }

        if (DotVec3(separation, separation) <= 1.0) {
            const double pi = std::acos(-1.0);
            return (1.0 + std::cos(pi * signed_distance / band_width)) / (2.0 * band_width);
        }
    }

    return 0.0;
}

void MPMMPMMonolithicFSI::LumpedLagrangeMultiplier() {

    std::vector<std::array<std::array<double, 3>, 3>> surface = this->BuildSolidInterface(this->solid_.nphi);
    const std::vector<std::array<std::array<double, 3>, 2>> fluid_domains = this->BuildFluidDomains();
    const std::map<std::array<int, 3>, std::vector<std::size_t>> fluid_bins = this->IndexFluidDomains(fluid_domains);
    interface_geometry::InterfaceSDF interface(std::move(surface), {dxy[0], dxy[1], dxy[2]});

    std::array<std::array<double, 3>, 6> quadrature_offset;
    GaussianDistribution(quadrature_offset);

    const double quadrature_volume = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nlm_lump_local);
    for (int m = 0; m < nelem; m++) {
        const std::array<int, 3> ijk = IndexToIJK(m, xyelem);

        std::array<double, 3> cell_center, xyg;
        for (int d = 0; d < 3; d++) { cell_center[d] = xymin[d] + dxy[d] * (double(ijk[d]) + 0.5); }

        for (int iz = 0; iz < npxye[2]; iz++) {
            xyg[2] = cell_center[2] + quadrature_offset[iz][2];
            for (int iy = 0; iy < npxye[1]; iy++) {
                xyg[1] = cell_center[1] + quadrature_offset[iy][1];
                for (int ix = 0; ix < npxye[0]; ix++) {
                    xyg[0] = cell_center[0] + quadrature_offset[ix][0];
                    const double area_density = this->WetInterfaceAreaDensity(xyg, interface, fluid_domains, fluid_bins);
                    if (area_density == 0.0) continue;

                    MakeSF(m, xyg, idimc, xynodec, ncm, nenode, sf, dsf);

                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        this->nlm_lump_local[nid] += sfi * quadrature_volume * area_density;
                    }
                }
            }
        }
    }

    this->nlm_lump = this->nlm_lump_local;

    NodeVarComm(this->nlm_lump, 0);

    return;
}

void MPMMPMMonolithicFSI::AssembleInterfaceSystem() {

    for (int n = 0; n < nodec; n++) {
        int ncol = 0;
        int ida = this->fsi_sys.FindIndex(n, n, ncol);
        const double lm = this->nlm_lump_local[n];
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[31]] -= lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[33]] -= lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[35]] -= lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[32]] += lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[34]] += lm;
        this->fsi_sys.amat[ida + this->fsi_sys.block_id[36]] += lm;
    }

    for (int n = 0; n < nodec; n++) {
        this->fsi_sys.b_rhs[n + nodec * 7] = this->nlm_lump[n] * (this->fluid_.ndispl[n + nuc] - this->solid_.ndispl[n + nuc]);
        this->fsi_sys.b_rhs[n + nodec * 8] = this->nlm_lump[n] * (this->fluid_.ndispl[n + nvc] - this->solid_.ndispl[n + nvc]);
        this->fsi_sys.b_rhs[n + nodec * 9] = this->nlm_lump[n] * (this->fluid_.ndispl[n + nwc] - this->solid_.ndispl[n + nwc]);
    }

    return;
}

void MPMMPMMonolithicFSI::AddLagrangeMultiplierToRHS(double af_coeff, const std::vector<int> &offsets) {

    for (int n = 0; n < nodec; n++) {
        this->fsi_sys.b_rhs[n + offsets[0]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + nuc];
        this->fsi_sys.b_rhs[n + offsets[1]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + nvc];
        this->fsi_sys.b_rhs[n + offsets[2]] += af_coeff * this->nlm_lump[n] * this->nlambda[n + nwc];
    }

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
        LocateLocalElement(xyp, cell);
        MakeSF(cell, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
        for (int ni = 0; ni < nenode; ni++) {
            const int nid = ncm[ni];
            const double sfi = sf[ni];
            this->solid_.StandardVarP2G(pid, nid, sfi, this->solid_.vol, final_phi);
        }
    }

    NodeVarComm(final_phi, 0);

    for (int n = 0; n < nodec; n++) { final_phi[n] = std::clamp(final_phi[n] / nvol[n], 0.0e0, 1.0e0); }

    interface_geometry::InterfaceSDF interface(this->BuildSolidInterface(final_phi), {dxy[0], dxy[1], dxy[2]});

    int corrected = 0;
    double maximum_shift = 0.0;
    for (int pid = 0; pid < this->fluid_.num; pid++) {
        const std::array<double, 3> original = this->fluid_.coord[pid];
        std::array<double, 3> surface_point, normal;
        double signed_distance;
        if (!interface.FindSurfacePoint(original, surface_point, signed_distance, normal) || signed_distance >= 0.0) continue;
        this->fluid_.coord[pid] = surface_point;
        maximum_shift = std::max(maximum_shift, -signed_distance);
        corrected++;
    }

    return;
}
