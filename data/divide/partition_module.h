#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include "../../module/material_point.h"
#include "../../module/bc.h"

using namespace std;

inline vector<int> nexyr1(3), nexyr2(3);
inline vector<int> nr1(3), nr2(3);
inline vector<int> nxyr1(3), nxyr2(3);
inline vector<double> drxy1(3), drxy2(3), bound12(3);
inline vector<int> inxyminc(3), inxymaxc(3), inxyminl(3), inxymaxl(3);

inline BoundaryCondition uwbcr, vwbcr, wwbcr, hpbcr,
                          usbcr, vsbcr, wsbcr,
                          uifbcr, vifbcr, wifbcr;

inline MaterialPoint wpr, spr;

void input_mesh_data (ifstream& infile);

void output_mesh_data (ofstream& outfile, int ir);

void partition_initial_dataset ();

void bc_renumber (BoundaryCondition& bcr, const BoundaryCondition& bc, vector<int> nxy_l,
                  vector<int> nxy_g, vector<int> nxymin, vector<int> nxymax, bool compute_fbcr = true);  

void point_renumber (int& local_num, vector<int>& global_id, MaterialPoint& point, int ir); 

void mesh_partition (int ir);

#endif