#include <cmath>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "material_point.h"
#include "./solid/solid_material_point.h"
#include "mesh.h"
#include "dataset.h"
#include "mpi_data.h"
#include "shape_function.h"
#include "contact.h"

using namespace std;

// void search_closest_point2node(int nid, MaterialPoint& point, double& distance) {

//     double near_point = -1;

//     double min_dist = 3.0e0 * dxy[0];
//     double radius = min_dist;
//     for (int pid = 0 ; pid < point.num; pid++) {
//         double dx = point.coord[0][pid] - xyc[0][nid];
//         double dy = point.coord[1][pid] - xyc[1][nid];
//         double dz = point.coord[2][pid] - xyc[2][nid];
//         double dist = sqrt(dx * dx + dy * dy + dz * dz);
//         if (dist > radius) {
//             continue;   // --- Skip if the distance is larger than the radius ---
//         }
//         if (dist < min_dist) {
//             min_dist = dist;
//             near_point = pid;
//         }
//     }

//     if (near_point == -1) {
//         distance = 0.0e0;  // --- Set distance to zero if no point is found ---
//         return;
//     }
    
//     double dx = (point.coord[0][near_point] - xyc[0][nid]) * point.nnormal[nid+nuc];
//     double dy = (point.coord[1][near_point] - xyc[1][nid]) * point.nnormal[nid+nvc];
//     double dz = (point.coord[2][near_point] - xyc[2][nid]) * point.nnormal[nid+nwc];
//     distance = dx + dy + dz;

//     return;
// }

// void search_fluid_solid_intf () {

//     // ---- Calculate the solid point normal vector at each node ----
//     cal_point_unit_normal(sp);
//     // ---- Calculate the water point normal vector at each node ----
//     VectorAssign(nodec*3, wp.nnormal);
//     for (int n = 0; n < nodec ; n++) {
//         wp.nnormal[n+nuc] = -sp.nnormal[n+nuc];
//         wp.nnormal[n+nvc] = -sp.nnormal[n+nvc];
//         wp.nnormal[n+nwc] = -sp.nnormal[n+nwc];
//     }
    
//     double dist_s, dist_w;
//     vector<double> cont_flag(nodec, 0.0e0);
//     for (int n = 0; n < nodec; n++) {
//         if (sp.nvof[n] > mtol && wp.nvof[n] > mtol) {
//             search_closest_point2node(n, sp, dist_s);
//             search_closest_point2node(n, wp, dist_w);
//             double distance = -(dist_s + dist_w);
//             if (distance < dxy[0] && distance > 0.0e0) {
//                 cont_flag[n] = 1.0e0 * dbc[n];
//             }
//         }
//     }

//     node_var_comm(cont_flag, 0);

//     fsi_intf.ibc = 0;
//     VectorAssign(nodec, fsi_intf.nbc);
//     for (int n = 0; n < nodec; n++) {
//         if (cont_flag[n] > 0.0e0) {
//             fsi_intf.nbc[fsi_intf.ibc++] = n;
//         }
//     }

//     // cout << fsi_intf.ibc << "\n";
//     // for (int i = 0; i < fsi_intf.ibc; i++) {
//     //    int nid = fsi_intf.nbc[i];
//     //    cout << "Contact node: " << setw(10) << i << setw(10) << nid << "\n"; 
//     // }

//     return;
// }