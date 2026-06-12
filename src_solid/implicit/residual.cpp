#include <cmath>
#include <vector>
#include <optional>
#include <string>
#include "../../module/bc.h"
#include "../../module/mesh.h"
#include "../../module/material_point.h"

void BCNRSet() {

    BCSetDt(nuc, sp.ndispl, usbc);
    BCSetDt(nvc, sp.ndispl, vsbc);
    BCSetDt(nwc, sp.ndispl, wsbc);

    return;
}