#include <cmath>
#include <vector>
#include <optional>
#include <string>
#include "../module/bc.h"
#include "../module/mesh.h"
#include "../module/solid/solid_material_point.h"
#include "../module/solid/implicit/implicit_mpm_solid.h"
#include "../module/fluid/FEM/stabilized_fem.h"
#include "block_fsi.h"

using namespace implicitmpm;
using namespace stabilizedfem;

// --- BCResidualSet & BCSet is used for self-defined solver ---
// --- BCPetscBCList is used for petsc solver ---

void FSIFluid::BCSet() {

    StabilizedFEM::BCSet();

    int num = this->fsi_.fsi_intf.ibc;
    if (num == 0) return; 
    for (int i = 0; i < num; i++) {
        int nid = this->fsi_.fsi_intf.nbc[i];
        this->nvel[nid+nuc] = this->fsi_.fsi_intf.fbc[i+0*num];
        this->nvel[nid+nvc] = this->fsi_.fsi_intf.fbc[i+1*num];
        this->nvel[nid+nwc] = this->fsi_.fsi_intf.fbc[i+2*num];
    }

    return;
}

void FSIFluid::BCResidualSet(std::vector<double>& rr) {

    StabilizedFEM::BCResidualSet(rr);

    this->fsi_.fsi_intf.BCSetZero(nuc, rr);
    this->fsi_.fsi_intf.BCSetZero(nvc, rr);
    this->fsi_.fsi_intf.BCSetZero(nwc, rr);

    return;
}

void FSIFluid::BuildPetscBCList(CrsMat& mat) {

    StabilizedFEM::BuildPetscBCList(mat);

    mat.AddBCComponent(this->fsi_.fsi_intf, nuc);
    mat.AddBCComponent(this->fsi_.fsi_intf, nvc);
    mat.AddBCComponent(this->fsi_.fsi_intf, nwc);

    return;
}
