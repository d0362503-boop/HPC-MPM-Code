#include <vector>
#include <cmath>
#include "../partition_module.h"
#include "../../../module/material_point.h"
#include "../../../module/bc.h"
#include "../../../module/mesh.h"

void partition_process (int ir) {

    // ---- CHANGE POINT RENUMBER DATA HERE (IF YOU NEED) ----
    point_renumber(spr.num, spr.id, sp, ir);
    // --------------------------------

    mesh_partition(ir);
 
    // ---- CHANGE BC RENUMBER DATA HERE (IF YOU NEED) ----
    bc_renumber(usbcr, usbc, xynodec, xynodecw, inxyminc, inxymaxc);
    bc_renumber(vsbcr, vsbc, xynodec, xynodecw, inxyminc, inxymaxc);
    bc_renumber(wsbcr, wsbc, xynodec, xynodecw, inxyminc, inxymaxc);

    bc_renumber(uwbcr, uwbc, xynodec, xynodecw, inxyminc, inxymaxc);
    bc_renumber(vwbcr, vwbc, xynodec, xynodecw, inxyminc, inxymaxc);
    bc_renumber(wwbcr, wwbc, xynodec, xynodecw, inxyminc, inxymaxc);
    bc_renumber(hpbcr, hpbc, xynodec, xynodecw, inxyminc, inxymaxc);

    // --- Inflow boundary condition ---
    bc_renumber(uifbcr, uifbc, xyelem, xyelemw, aelemmin, aelemmax, false);
    bc_renumber(vifbcr, vifbc, xyelem, xyelemw, aelemmin, aelemmax, false);
    bc_renumber(wifbcr, wifbc, xyelem, xyelemw, aelemmin, aelemmax, false);
    // --------------------------------

    return;
}
