#include <vector>
#include <cmath>
#include "../../module/material_point.h"
#include "../../module/mpi_data.h"
#include "../../module/mesh.h"
#include "../../module/bc.h"
#include "../../module/dataset.h"
#include "../../module/data_IO.h"
#include "partition_module.h"

using namespace std;

void input_mesh_data (ifstream& infile) {

    for (int i = 0; i < 3; i++) infile >> idimc[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> xyminw[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> xymaxw[i];
    infile.ignore(1000, '\n');

    infile >> nelem; for (int i = 0; i < 3; i++) infile >> xyelemw[i];
    infile.ignore(1000, '\n');

    infile >> node;  for (int i = 0; i < 3; i++) infile >> xynodew[i];
    infile.ignore(1000, '\n');

    infile >> nodec; for (int i = 0; i < 3; i++) infile >> xynodecw[i];
    infile.ignore(1000, '\n');

    return;
}

void output_mesh_data (ofstream& outfile, int ir) {

    outfile << setw(15) << ir << "     --- myrank ---" << "\n";

    outfile << setw(15) << nelem  << 
               setw(15) << xyelem[0] << 
               setw(15) << xyelem[1] << 
               setw(15) << xyelem[2] << "     --- rank element ---" << "\n";

    outfile << setw(15) << node  << 
               setw(15) << xynode[0] << 
               setw(15) << xynode[1] << 
               setw(15) << xynode[2] << "     --- rank node ---" << "\n";

    outfile << setw(15) << nodec  << 
               setw(15) << xynodec[0] << 
               setw(15) << xynodec[1] << 
               setw(15) << xynodec[2] << "     --- rank control point ---" << "\n";

    for (int i = 0; i < 3; i++) outfile << setw(15) << nxyr[i];
    outfile << "     --- rank number ---" << "\n";

    for (int i = 0; i < 3; i++) outfile << setw(15) << xymin[i];
    outfile << "     --- rank minimum coordinate ---" << "\n";

    for (int i = 0; i < 3; i++) outfile << setw(15) << xymax[i];
    outfile << "     --- rank maximum coordinate ---" << "\n";

    for (int i = 0; i < 3; i++) outfile << setw(15) << aelemmin[i];
    outfile << "     --- rank minimum element ID ---" << "\n";

    for (int i = 0; i < 3; i++) outfile << setw(15) << aelemmax[i];
    outfile << "     --- rank maximum element ID ---" << "\n";

    outfile << setw(15) << isb << "     --- overlap rank number ---" << "\n";
    if (isb != 0) {
        output_oned_vector(outfile, isb, naid);

        outfile << setw(15) << isubc << "     --- overlap control point number ---" << "\n";
        output_oned_vector(outfile, isb, nsbc);

        output_oned_vector(outfile, isubc, nsubc);

        output_oned_vector(outfile, nodec, dbc);

        outfile << setw(15) << isubl << "     --- overlap node number ---" << "\n";
        output_oned_vector(outfile, isb, nsbl);

        output_oned_vector(outfile, isubl, nsubl);

        output_oned_vector(outfile, node, dbl);
    }
    return;
}

void partition_initial_dataset () {

    for (int i = 0; i < 3; i++) {
        nexyr1[i] = xyelemw[i] / nxyr[i];
        nexyr2[i] = nexyr1[i] + 1;
        nr2[i] = xyelemw[i] % nxyr[i];
        nr1[i] = nxyr[i] - nr2[i];
        dxy[i] = (xymaxw[i] - xyminw[i]) / double(xyelemw[i]);
        drxy1[i] = (xymaxw[i] - xyminw[i] - dxy[i] * double(nr2[i])) / double(nxyr[i]);
        drxy2[i] = drxy1[i] + dxy[i];
        nxyr1[i] = nexyr1[i] + 1;
        nxyr2[i] = nexyr2[i] + 1;
        bound12[i] = xyminw[i] + drxy1[i] * double(nr1[i]);
    }
    return;
}

void bc_renumber (BoundaryCondition& bcr, const BoundaryCondition& bc, vector<int> nxy_l,
                  vector<int> nxy_g, vector<int> nxymin, vector<int> nxymax, bool compute_fbcr) {               

    bcr.nbc.resize(bc.ibc);
    if (compute_fbcr == true) {
        bcr.fbc.resize(bc.ibc);
    }

    vector<int> inxy(3), inxyr(3);                
    int nx = nxy_g[0];
    int ny = nxy_g[1];
    int nz = nxy_g[2];

    bcr.ibc = 0;
    for (int i = 0; i < bc.ibc; ++i) {
        int n = bc.nbc[i];
        inxy[2] = n / (nx * ny);
        inxy[1] = (n - inxy[2] * (nx * ny)) / nx;
        inxy[0] = n - inxy[2] * (nx * ny) - inxy[1] * nx;

        if (inxy[0] >= nxymin[0] && inxy[0] <= nxymax[0] &&
            inxy[1] >= nxymin[1] && inxy[1] <= nxymax[1] &&
            inxy[2] >= nxymin[2] && inxy[2] <= nxymax[2]) {
            for (int j = 0; j < 3; ++j) {
                inxyr[j] = inxy[j] - nxymin[j];
            }
            int inr = nxy_l[0] * nxy_l[1] * inxyr[2] + nxy_l[0] * inxyr[1] + inxyr[0];
            bcr.nbc[bcr.ibc] = inr;
            if (compute_fbcr == true) {
                bcr.fbc[bcr.ibc] = bc.fbc[i]; 
            }
            ++bcr.ibc;
        }
    }

    return;
}

void point_renumber (int& local_num, vector<int>& global_id, MaterialPoint& point, int ir) { 

     vector<int> irxy(3);
     global_id.resize(point.num);
     int nrx = nxyr[0], nry = nxyr[1], nrz = nxyr[2];
     
     local_num = 0;
     for (int ip = 0; ip < point.num; ip++) {
        for (int i = 0; i < 3; i++) {
            if (point.coord[ip][i] < bound12[i]) {
                irxy[i] = floor((point.coord[ip][i] - xyminw[i]) / drxy1[i]);
            }
            else {
                irxy[i] = nr1[i] + floor((point.coord[ip][i] - bound12[i]) / drxy2[i]);
            }
        }
        int irp = nrx * nry * irxy[2] + nrx * irxy[1] + irxy[0];
        if (irp == ir) {
            global_id[local_num] = ip;
            local_num++;
        } 
    }

    return;
}

void mesh_partition (int ir) {

    int nrx = nxyr[0], nry = nxyr[1], nrz = nxyr[2];
    int nex = xyelemw[0], ney = xyelemw[1], nez = xyelemw[2];
    int nxl = xynodew[0], nyl = xynodew[1], nzl = xynodew[2];
    int nxc = xynodecw[0], nyc = xynodecw[1], nzc = xynodecw[2];

    vector<int> irxy(3);
    irxy[2] = ir / (nrx * nry);
    irxy[1] = (ir - irxy[2] * (nrx * nry)) / nrx;
    irxy[0] = ir - irxy[2] * (nrx * nry) - irxy[1] * nrx;

    isb = 1;

    vector<vector<int>> ijkl(6,vector<int>(3)), ijkc(6,vector<int>(3));
    for (int ix = 0; ix < 3; ix++) {

        if (irxy[ix] < nr1[ix]) {
            xymin[ix] = xyminw[ix] + drxy1[ix] * double(irxy[ix]);
            xymax[ix] = xymin[ix] + drxy1[ix];
            aelemmin[ix] = nexyr1[ix] * irxy[ix];
            aelemmax[ix] = aelemmin[ix] + nexyr1[ix] - 1;
            xyelem[ix] = nexyr1[ix];
        }  
        else {
            xymin[ix] = bound12[ix] + drxy2[ix] * double(irxy[ix] - nr1[ix]);
            xymax[ix] = xymin[ix] + drxy2[ix];
            aelemmin[ix] = nexyr1[ix] * nr1[ix] + nexyr2[ix] * (irxy[ix] - nr1[ix]);
            aelemmax[ix] = aelemmin[ix] + nexyr2[ix] - 1;
            xyelem[ix] = nexyr2[ix];
        }

        xynode[ix] = xyelem[ix] + 1;
        xynodec[ix] = xyelem[ix] + idimc[ix];

        if (nxyr[ix] == 1) {
            int temp_ijkl[] = {-1, -1, 0, xynode[ix]-1, -1, -1};
            int temp_ijkc[] = {-1, -1, 0, xynodec[ix]-1, -1, -1};
            for (int i = 0; i < 6; i++) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
        } 
        else if (irxy[ix] == 0) {  // Note: Adjusted for 0-based indexing
            int temp_ijkl[] = {-1, -1, 0, xynode[ix]-1, xynode[ix]-1, xynode[ix]-1};
            int temp_ijkc[] = {-1, -1, 0, xynodec[ix]-1, xynodec[ix]-idimc[ix], xynodec[ix]-1};
            for (int i = 0; i < 6; i++) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            isb *= 2;
        }
        else if (irxy[ix] == nxyr[ix] - 1) {  // Adjusted for 0-based indexing
            int temp_ijkl[] = {0, 0, 0, xynode[ix]-1, -1, -1};
            int temp_ijkc[] = {0, idimc[ix]-1, 0, xynodec[ix]-1, -1, -1};
            for (int i = 0; i < 6; i++) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            isb *= 2;
        }
        else {
            int temp_ijkl[] = {0, 0, 0, xynode[ix]-1, xynode[ix]-1, xynode[ix]-1};
            int temp_ijkc[] = {0, idimc[ix]-1, 0, xynodec[ix]-1, xynodec[ix]-idimc[ix], xynodec[ix]-1};
            for (int i = 0; i < 6; i++) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            isb *= 3;
        }

        inxyminl[ix] = aelemmin[ix];
        inxymaxl[ix] = aelemmax[ix] + 1;
        inxyminc[ix] = aelemmin[ix];
        inxymaxc[ix] = aelemmax[ix] + idimc[ix];
    }

    node = xynode[0] * xynode[1] * xynode[2];
    nodec = xynodec[0] * xynodec[1] * xynodec[2];
    nelem = xyelem[0] * xyelem[1] * xyelem[2];

    isb--; 
    naid.resize(isb);
    nsbc.resize(isb);
    nsbl.resize(isb);
    nsubc.resize(nodec*isb); 
    nsubl.resize(node*isb);
    VectorAssign(nodec, dbc, 1.0e0);
    VectorAssign(node, dbl, 1.0e0);

    isubc = 0;
    isubl = 0;
    isb = 0;
    for (int kk = 0; kk < 3; kk++) {
        int ksc = ijkc[2 * kk][2]; 
        int kec = ijkc[2 * kk + 1][2];
        int ksl = ijkl[2 * kk][2]; 
        int kel = ijkl[2 * kk + 1][2];
        if (ksc == -1) continue;
        for (int jj = 0; jj < 3; jj++) {
            int jsc = ijkc[2 * jj][1]; 
            int jec = ijkc[2 * jj + 1][1];
            int jsl = ijkl[2 * jj][1]; 
            int jel = ijkl[2 * jj + 1][1];
            if (jsc == -1) continue;
            for (int ii = 0; ii < 3; ii++) {
                if (ii * jj * kk == 1) continue;  //Skip rank own 
                int isc = ijkc[2 * ii][0]; 
                int iec = ijkc[2 * ii + 1][0];
                int isl = ijkl[2 * ii][0]; 
                int iel = ijkl[2 * ii + 1][0];
                if (isc == -1) continue;
                naid[isb] = ir + (ii - 1) + (jj - 1) * nxyr[0] + (kk - 1) * nxyr[0] * nxyr[1]; 
                nsbc[isb] = (iec - isc + 1) * (jec - jsc + 1) * (kec - ksc + 1);
                nsbl[isb] = (iel - isl + 1) * (jel - jsl + 1) * (kel - ksl + 1);
                for (int k = ksl; k <= kel; k++) {
                    for (int j = jsl; j <= jel; j++) {
                        for (int i = isl; i <= iel; i++) {
                          int nn = i + j * xynode[0] + k * xynode[0] * xynode[1];
                          nsubl[isubl++] = nn;
                          dbl[nn] += 1.0e0;
                        }
                    }
                }
                for (int k = ksc; k <= kec; k++) {
                    for (int j = jsc; j <= jec; j++) {
                        for (int i = isc; i <= iec; i++) {
                          int nn = i + j * xynodec[0] + k * xynodec[0] * xynodec[1];
                          nsubc[isubc++] = nn;
                          dbc[nn] += 1.0e0;
                        }
                    }
                }
                isb++;
            }
        }
    }

    return;
}