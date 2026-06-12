#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include "dataset.h"
#include "mpi_data.h"
#include "material_point.h"
#include "shape_function.h"


// std::vector<std::vector<double>> particle_shifting_SPH_like(MaterialPoint& point, Particle_Movement& point_move) {
    
//     VectorAssign(3, point.num, spring_force);

//     // --- adjacent ghost point communication ---
//     std::vector<int> mxy(3), checkmp;
//     checkmp.assign(point.num, -1);
//     point_move.nsp.assign(isb, 0);
//     point_move.nrp.assign(isb, 0);
//     for (int ip = 0; ip < point.num; ip++) {
//         double radius = cbrt(0.75e0 * point.vol[ip] / M_PI);
//         for (int i = 0; i < 3; i++) {
//             if (abs(point.coord[i][ip]-xymin[i]) < radius) {
//                 mxy[i] = -1;
//             } else if (abs(point.coord[i][ip]-xymax[i]) < radius) {
//                 mxy[i] = 1;
//             }
//         }

//         int idsb = myrank + mxy[0] + mxy[1] * nxyr[0] 
//                           + mxy[2] * nxyr[0] * nxyr[1];
        
//         for (int i = 0; i < isb; i++) {
//             if (naid[i] == idsb) {
//                 point_move.nsp[i]++;
//                 checkmp[ip] = i;
//                 break;
//             }
//         }
//     }

//     std::vector<int> idmploc(isb);
//     point_move.nmps = 0;
//     for (int i = 0; i < isb; i++) {
//         idmploc[i] = point_move.nmps;
//         point_move.nmps += point_move.nsp[i];
//     }

//     VectorAssign(point_move.nmps, point_move.idmp);

//     for (int ip = 0; ip < point.num; ip++) {
//         int i = checkmp[ip];
//         if (i != -1) { // Move
//             point_move.idmp[idmploc[i]++] = ip;
//         }
//     }

//     std::vector<MPI_Request> irqs(isb);
//     std::vector<MPI_Request> irqr(isb);  
//     MPI_Status status;
//     for (int i = 0; i < isb; i++) {
//         int ncomid = naid[i];
//         MPI_Isend(&point_move.nsp[i], 1, MPI_INT, ncomid, 1, MPI_COMM_WORLD, &irqs[i]);
//         MPI_Irecv(&point_move.nrp[i], 1, MPI_INT, ncomid, 1, MPI_COMM_WORLD, &irqr[i]);
//     }
//     for (int i = 0; i < isb; i++) {
//         MPI_Wait(&irqs[i], &status);
//         MPI_Wait(&irqr[i], &status);
//     }

//     point_move.nrps = 0;
//     for (int i = 0; i < isb; i++) {
//         point_move.nrps += point_move.nrp[i];
//     }

//     MaterialPoint ghost_point;

//     ghost_point.num = point_move.nrps;
//     VectorAssign(3, ghost_point.num, ghost_point.coord);
    
//     int sip = 0, idm = 3; 
//     std::vector<double> bufs(point_move.nmps*idm), bufr(point_move.nrps*idm);
//     for (int i = 0; i < isb; i++) {
//         for (int ip = 0; ip < point_move.nsp[i]; ip++) {
//             int ipsip = ip + sip;
//             int ipsip2 = idm * ipsip;
//             for (int j = 0; j < idm; j++) {
//                 bufs[j+ipsip2] = point.coord[j][point_move.idmp[ipsip]];
//             }
//         }
//         sip += point_move.nsp[i];
//     }
    
//     point_var_sendrecv(bufs, bufr, point_move, idm);

//     for (int j = 0; j < point_move.nrps; j++) {
//         for (int i = 0; i < idm; i++) {
//             ghost_point.coord[i][j] = bufr[i+j*idm];
//         }
//     }
//     // ------

//     std::vector<int> num_near_point;
//     std::vector<std::vector<double>> near_point_coord;
//     VectorAssign(point.num, num_near_point);
//     VectorAssign(3, point.num+point_move.nrps, near_point_coord);

//     for (int ip = 0; ip < point.num; ip++){
//         double radius = cbrt(0.75e0 * point.vol[ip] / M_PI);
//         for (int jp = 0; jp < point.num; jp++) {
//             double dx = point.coord[0][jp] - point.coord[0][ip];
//             double dy = point.coord[1][jp] - point.coord[1][ip];
//             double dz = point.coord[2][jp] - point.coord[2][ip];
//             double norm = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
//             if (jp != ip && (norm < radius && norm > 0)) {
//                 int n = num_near_point[ip];
//                 near_point_coord[0][n] = point.coord[0][jp];
//                 near_point_coord[1][n] = point.coord[1][jp];
//                 near_point_coord[2][n] = point.coord[2][jp];
//                 num_near_point[ip]++;
//             } 
//         }
//         for (int jp = 0; jp < ghost_point.num; jp++) {
//             double dx = ghost_point.coord[0][jp] - point.coord[0][ip];
//             double dy = ghost_point.coord[1][jp] - point.coord[1][ip];
//             double dz = ghost_point.coord[2][jp] - point.coord[2][ip];
//             double norm = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
//             if (norm < radius && norm > 0) {
//                 int n = num_near_point[ip];
//                 near_point_coord[0][n] = ghost_point.coord[0][jp];
//                 near_point_coord[1][n] = ghost_point.coord[1][jp];
//                 near_point_coord[2][n] = ghost_point.coord[2][jp];
//                 num_near_point[ip]++;
//             } 
//         }

//         const double a = 50.0e0;
//         for (int jp = 0; jp < num_near_point[ip]; jp++) {
//             double dx = near_point_coord[0][jp] - point.coord[0][ip];
//             double dy = near_point_coord[1][jp] - point.coord[1][ip];
//             double dz = near_point_coord[2][jp] - point.coord[2][ip];
//             double norm = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
//             double smooth_weight = (norm < mtol) ? 0.0e0 : (1.0e0 - pow(norm, 2) / pow(radius, 2));
//             spring_force[0][ip] += dx / norm * smooth_weight;
//             spring_force[1][ip] += dy / norm * smooth_weight;
//             spring_force[2][ip] += dz / norm * smooth_weight;
//         }
//         for (int i = 0; i < 3; i++) {
//             spring_force[i][ip] *= -a * radius * dt; 
//         }
//     }

//     return;
// }
