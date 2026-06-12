#include <cmath>
#include <vector>
#include "../../bc.h"
#include "../../mesh.h"
#include "../../material_point.h"
#include "../../dataset.h"

void SolveSolid() {

    VectorAssign(nodec*3, sp.naccel);    
    for (int n = 0; n < nodec; n++) {
        if (sp.nmass[n] > mtol) {
            sp.naccel[n+nuc] = sp.nforce[n+nuc] / sp.nmass[n];
            sp.naccel[n+nvc] = sp.nforce[n+nvc] / sp.nmass[n];
            sp.naccel[n+nwc] = sp.nforce[n+nwc] / sp.nmass[n];
        }
        else {
            sp.naccel[n+nuc] = 0.0e0;
            sp.naccel[n+nvc] = 0.0e0;
            sp.naccel[n+nwc] = 0.0e0;
        }
    }

    BCSetZero(nuc, sp.naccel, usbc);
    BCSetZero(nvc, sp.naccel, vsbc);
    BCSetZero(nwc, sp.naccel, wsbc);
    BCSetZero(nuc, sp.naccel, rigid_bc);
    BCSetZero(nvc, sp.naccel, rigid_bc);
    BCSetZero(nwc, sp.naccel, rigid_bc);

    return;
}