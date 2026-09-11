#include "module/bc.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"
#include "module/solver/crsmat.h"
#include "work/src_fsi/MPM_MPM/block_fsi.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

using namespace mpmmpmblockfsi;

void FSISolid::AssembleSystem(const std::vector<double> &naccel_k, //
                              const std::vector<double> &nvel_k,   //
                              std::vector<std::array<double, 6>> &stress_k) {

    implicitmpm::ImplicitSolidMPM::AssembleSystem(naccel_k, nvel_k, stress_k);

    for (int n = 0; n < nodec; n++) {
        this->SM_.b_rhs[n + nuc] += this->lm_lumped[n] * this->fsi_.nfsi_force[n + nuc];
        this->SM_.b_rhs[n + nvc] += this->lm_lumped[n] * this->fsi_.nfsi_force[n + nvc];
        this->SM_.b_rhs[n + nwc] += this->lm_lumped[n] * this->fsi_.nfsi_force[n + nwc];
    }

    return;
}
