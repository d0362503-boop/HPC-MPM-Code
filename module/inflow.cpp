#include <algorithm>
#include <mpi.h>

#include "module/dataset.h"
#include "module/material_point.h"

bool MaterialPoint::InflowMeshisFilled(const BoundaryCondition &infbc) {
    const RunParams &params = this->params_;

    const int target_per_cell = params.npxye[0] * params.npxye[1] * params.npxye[2];

    bool isfilled = false;
    static double fill_ratio_tol = 0.95e0;

    for (int m = 0; m < infbc.ibc; m++) {
        int ie = infbc.nbc[m];
        double fill_ratio = this->numep[ie] / static_cast<double>(target_per_cell);
        if (fill_ratio > fill_ratio_tol) isfilled = true;
    }

    return isfilled;
}

void MaterialPoint::GenerateInflowParticles(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) {

    if (!this->infbc_isfilled) { this->infbc_isfilled = InflowMeshisFilled(infbc); }

    if (this->infbc_isfilled) {
        this->GenerateInflowParticlesFilledMesh(dir, ifp, infbc);
    } else {
        this->GenerateInflowParticlesEmptyMesh(dir, ifp, infbc);
    }

    return;
}

void MaterialPoint::AssignUniqueInflowIds(MaterialPoint &ifp) const {

    int local_max_id = -1;
    if (!this->id.empty()) { local_max_id = *std::max_element(this->id.begin(), this->id.end()); }

    int global_max_id = -1;
    MPI_Allreduce(&local_max_id, &global_max_id, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    int local_new_count = ifp.num;
    int id_offset = 0;
    MPI_Exscan(&local_new_count, &id_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (myrank == 0) { id_offset = 0; }

    const int next_id = global_max_id + 1 + id_offset;
    for (int i = 0; i < ifp.num; ++i) { ifp.id[i] = next_id + i; }

    return;
}
