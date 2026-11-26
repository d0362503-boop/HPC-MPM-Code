#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solid/solid_material_point.h"
#include "work/src_fsi/MPM_FEM/block_fsi.h"
#include <vector>

using namespace implicitmpm;
using namespace stabilizedfem;
using namespace mpmfemblockfsi;

// --- BCResidualSet & BCSet is used for self-defined solver ---
// --- BCPetscBCList is used for petsc solver ---

void FSIFluid::BCNRSet(const TimeState &ts) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    StabilizedFEM::BCNRSet(ts);

    int num = this->fsi_.fsi_intf.ibc;
    if (num == 0) return;
    for (int i = 0; i < num; i++) {
        int nid = this->fsi_.fsi_intf.nbc[i];
        this->nvel_k[nid + cmesh.nuc] = this->fsi_.fsi_intf.fbc[i + 0 * num];
        this->nvel_k[nid + cmesh.nvc] = this->fsi_.fsi_intf.fbc[i + 1 * num];
        this->nvel_k[nid + cmesh.nwc] = this->fsi_.fsi_intf.fbc[i + 2 * num];
    }

    return;
}

void FSIFluid::BCResidualSet(std::vector<double> &rr) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    StabilizedFEM::BCResidualSet(rr);

    this->fsi_.fsi_intf.BCSetZero(cmesh.nuc, rr);
    this->fsi_.fsi_intf.BCSetZero(cmesh.nvc, rr);
    this->fsi_.fsi_intf.BCSetZero(cmesh.nwc, rr);

    return;
}

void FSIFluid::BuildPetscBCList(CrsMat &mat) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    StabilizedFEM::BuildPetscBCList(mat);

    mat.AddBCComponent(this->fsi_.fsi_intf, cmesh.nuc);
    mat.AddBCComponent(this->fsi_.fsi_intf, cmesh.nvc);
    mat.AddBCComponent(this->fsi_.fsi_intf, cmesh.nwc);

    return;
}
