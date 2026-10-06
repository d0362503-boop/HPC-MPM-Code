#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <string>
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

void MPMMPMMonolithicFSI::SolveFSISystem() {

    this->LumpedLagrangeMultiplier();

    this->BuildActiveDOFs();

    this->fluid_.MakeNSStabCoeff(this->fluid_.nvel);

    std::vector<std::array<double, 6>> stress_k = this->solid_.InitializeNRStress();

    std::vector<double> nvel_f(nodec * 3), nvel_s(nodec * 3);
    std::vector<double> naccel_f(nodec * 3), naccel_s(nodec * 3);
    std::array<double, 4> initial_norm{};

    VectorAssign(nodec * 3, this->fluid_.ndispl);
    VectorAssign(nodec * 3, this->solid_.ndispl);
    VectorAssign(nodec * 3, this->nlambda);
    VectorAssign(nodec, this->fluid_.npres);
    VectorAssign(nodec * 10, this->fsi_sys.x_lhs);

    for (int NR_it = 0; NR_it <= this->max_NR_it; NR_it++) {
        VectorAssign(this->fsi_sys.nmata, this->fsi_sys.amat);
        VectorAssign(nodec * 10, this->fsi_sys.b_rhs);

        this->fluid_.BCNRSet();
        this->fluid_.ComputeNodeVelAccelFromDispl(nvel_f, naccel_f);

        this->solid_.BCNRSet();
        this->solid_.ComputeNodeVelAccelFromDispl(nvel_s, naccel_s);

        this->AssembleFluidSystem(nvel_f, naccel_f);
        this->AssembleSolidSystem(nvel_s, naccel_s, stress_k);
        this->AssembleInterfaceSystem();

        int solver_it = this->SolveSystem(NR_it);

        this->UpdateNRIncrement();

        if (this->CheckNRConvergence(initial_norm, NR_it, solver_it)) { break; }
    }

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

void MPMMPMMonolithicFSI::LumpedLagrangeMultiplier() {

    std::vector<std::array<std::array<double, 3>, 3>> surface = this->BuildSolidInterface();
    std::vector<std::array<std::array<double, 3>, 2>> fluid_domains = this->BuildFluidDomains();
    interface_geometry::InterfaceSDF interface(std::move(surface), std::move(fluid_domains));

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
                    const double area_density = interface.AreaDensity(xyg);
                    if (area_density == 0.0) continue;

                    MakeSF(m, xyg, idimc, xynodec, ncm, nenode, sf, dsf);

                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        this->nlm_lump_local[nid] += sf[ni] * quadrature_volume * area_density;
                    }
                }
            }
        }
    }

    this->nlm_lump = this->nlm_lump_local;

    NodeVarComm(this->nlm_lump, 0);

    return;
}

void MPMMPMMonolithicFSI::OutputInterfaceBalance() const {
    std::array<double, 7> balance{};
    for (int n = 0; n < nodec; ++n) {
        const double area = this->nlm_lump_local[n];
        balance[0] += area;
        for (int d = 0; d < 3; ++d) {
            const double traction_force = area * this->nlambda[n + d * nodec];
            balance[1 + d] += this->fluid_.integrator_.alpha_f * traction_force;
            balance[4 + d] -= this->solid_.integrator_.alpha_f * traction_force;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, balance.data(), balance.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (myrank == 0) {
        std::ostringstream message;
        message << std::scientific << std::setprecision(9) << "Interface_balance: " << istep;
        for (double value : balance) { message << ' ' << value; }
        std::cout << message.str() << '\n';
    }
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

std::vector<std::array<std::array<double, 3>, 3>> MPMMPMMonolithicFSI::BuildSolidInterface() const {

    // Consistent tetrahedral subdivision across cell faces.
    const std::array<std::array<int, 4>, 6> tetrahedra{
        {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}}};
    std::vector<std::array<std::array<double, 3>, 3>> surface;
    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    const int nx = npxye[0] + 1, ny = npxye[1] + 1, nz = npxye[2] + 1;
    std::vector<std::array<double, 3>> vertices(nx * ny * nz);
    std::vector<double> level(vertices.size());

    for (int m = 0; m < nelem; m++) {
        double minimum_phi = 1.0, maximum_phi = 0.0;
        for (int nid : ncc[m]) {
            minimum_phi = std::min(minimum_phi, this->solid_.nphi[nid]);
            maximum_phi = std::max(maximum_phi, this->solid_.nphi[nid]);
        }
        if (minimum_phi >= 0.5 || maximum_phi <= 0.5) continue;

        const std::array<int, 3> ijk = IndexToIJK(m, xyelem);
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    const int vertex_id = i + nx * (j + ny * k);
                    const std::array<int, 3> offset{i, j, k};
                    std::array<double, 3> &xyg = vertices[vertex_id];
                    for (int d = 0; d < 3; d++) { xyg[d] = xymin[d] + dxy[d] * (ijk[d] + double(offset[d]) / npxye[d]); }

                    MakeSF(m, xyg, idimc, xynodec, ncm, nenode, sf, dsf);

                    level[vertex_id] = 0.5;
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        level[vertex_id] -= sf[ni] * this->solid_.nphi[nid];
                    }
                }
            }
        }

        for (int k = 0; k < npxye[2]; k++) {
            for (int j = 0; j < npxye[1]; j++) {
                for (int i = 0; i < npxye[0]; i++) {
                    std::array<int, 8> cube;
                    for (int corner = 0; corner < 8; corner++) {
                        cube[corner] = i + (corner & 1) + nx * (j + ((corner >> 1) & 1) + ny * (k + (corner >> 2)));
                    }

                    for (const std::array<int, 4> &tetrahedron : tetrahedra) {
                        std::array<std::array<double, 3>, 4> positions;
                        std::array<double, 4> values;
                        for (int corner = 0; corner < 4; corner++) {
                            const int vertex_id = cube[tetrahedron[corner]];
                            positions[corner] = vertices[vertex_id];
                            values[corner] = level[vertex_id];
                        }

                        this->AppendSolidContour(positions, values, surface);
                    }
                }
            }
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

void MPMMPMMonolithicFSI::UpdateNRIncrement() {

    for (int n = 0; n < nodec * 3; n++) {
        this->fluid_.ndispl[n] += this->fsi_sys.x_lhs[n];
        this->solid_.ndispl[n] += this->fsi_sys.x_lhs[n + nodec * 4];
        this->nlambda[n] += this->fsi_sys.x_lhs[n + nodec * 7];
    }

    for (int n = 0; n < nodec; n++) { this->fluid_.npres[n] += this->fsi_sys.x_lhs[n + npc]; }

    return;
}

bool MPMMPMMonolithicFSI::CheckNRConvergence(std::array<double, 4> &initial_norm, int NR_it, int solver_it) {

    const std::array<double, 4> absolute_tol = {1.0e-8, 1.0e-10, 1.0e-8, 1.0e-8};

    // Compute b - A*x in the scaled PETSc system.
    Vec linear_residual, local_residual;
    VecDuplicate(this->fsi_sys.petsc_b, &linear_residual);
    MatMult(this->fsi_sys.petsc_mat, this->fsi_sys.petsc_x, linear_residual);
    VecAYPX(linear_residual, -1.0, this->fsi_sys.petsc_b); // residual = b - A*x

    // Gather residuals in the same local ordering as x_lhs.
    VecDuplicate(this->fsi_sys.seq_x, &local_residual);
    VecScatterBegin(this->fsi_sys.scatter_to_all, linear_residual, local_residual, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(this->fsi_sys.scatter_to_all, linear_residual, local_residual, INSERT_VALUES, SCATTER_FORWARD);
    const PetscScalar *scaled_residual;
    VecGetArrayRead(local_residual, &scaled_residual);

    // stats[0]: linearized fluid residual [N].
    // stats[1]: linearized continuity residual [m^3/s].
    // stats[2]: linearized solid residual [N].
    // stats[3]: displacement increment difference [m].
    // stats[4..7]: corresponding active-component counts.
    std::array<double, 8> stats{};

    for (int n : this->fsi_sys.owned_natural_ids) {
        for (int d = 0; d < 10; d++) {
            const int index = n + d * nodec;
            if (this->fixed_dof[index]) continue;

            int field;
            if (d < 3) { // Fluid momentum
                field = 0;
            } else if (d == 3) { // Fluid continuity
                field = 1;
            } else if (d < 7) { // Solid momentum
                field = 2;
            } else { // Interface displacement
                field = 3;
            }

            double residual;
            if (d < 7) {
                // Undo row scaling; it equals column scaling for these seven fields.
                const double row_scale = this->column_scale[index];
                residual = scaled_residual[index] / row_scale;
            } else {
                const int displacement_index = n + (d - 7) * nodec;
                residual = this->solid_.ndispl[displacement_index] - this->fluid_.ndispl[displacement_index];
            }
            stats[field] += residual * residual;
            stats[field + 4] += 1.0;
        }
    }

    VecRestoreArrayRead(local_residual, &scaled_residual);
    VecDestroy(&local_residual);
    VecDestroy(&linear_residual);

    MPI_Allreduce(MPI_IN_PLACE, stats.data(), 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // All four RMS residuals must be finite and satisfy their thresholds.
    // Fields 0..2: max(absolute_tol, 1e-4 * NR0 post-solve residual RMS).
    // Field 3: displacement-increment RMS <= 1e-8 m, without interface weights.
    bool converged = true;
    for (int field = 0; field < 4; field++) {
        stats[field] = std::sqrt(stats[field] / std::max(1.0, stats[field + 4]));
        if (NR_it == 0 && field < 3) { initial_norm[field] = stats[field]; }
        const double tolerance = field == 3 ? absolute_tol[field] : std::max(absolute_tol[field], 1.0e-4 * initial_norm[field]);
        converged = converged && stats[field] <= tolerance;
    }

    if (converged) {
        if (myrank == 0) {
            std::cout << "Monolithic_converge: " << std::setw(15) << NR_it << std::setw(15) << solver_it //
                      << std::scientific << std::setw(15) << stats[0] << std::setw(15) << stats[1]       //
                      << std::setw(15) << stats[2] << std::setw(15) << stats[3] << "\n";
        }
        return true;
    }

    if (myrank == 0) {
        std::cout << "Monolithic_NR: " << std::setw(15) << NR_it << std::setw(15) << solver_it //
                  << std::scientific << std::setw(15) << stats[0] << std::setw(15) << stats[1] //
                  << std::setw(15) << stats[2] << std::setw(15) << stats[3] << "\n";
    }

    if (NR_it == this->max_NR_it) {
        if (myrank == 0) { std::cout << "Monolithic Newton did not converge" << "\n"; }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return false;
}
