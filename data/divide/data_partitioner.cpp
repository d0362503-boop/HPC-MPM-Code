#include "data_partitioner.h"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"

int DataPartitioner::Run() {
    std::cout << " ---- Start dividing " << this->CaseName() << " data ----"
              << "\n";

    std::ifstream para_infile = OpenInputFile("para_input_data.txt");
    this->LoadPartitionParameters(para_infile);

    std::ifstream grid_infile = OpenInputFile(this->grid_input_file_);
    this->InputMeshData(grid_infile);
    this->LoadBoundaryData(grid_infile);

    std::ifstream point_infile = OpenInputFile(this->point_input_file_);
    this->LoadPointData(point_infile);

    this->PartitionInitialDataset();
    this->CreateOutputDirectories();

    const int rank_count = nxyr[0] * nxyr[1] * nxyr[2];
    for (int rank_id = 0; rank_id < rank_count; ++rank_id) {
        this->PartitionProcess(rank_id);

        std::ofstream grid_outfile =
            OpenOutputFile(this->grid_output_prefix_ + std::to_string(rank_id) + ".txt");
        this->OutputMeshData(grid_outfile, rank_id);
        this->WriteBoundaryData(grid_outfile);

        std::ofstream point_outfile =
            OpenOutputFile(this->point_output_prefix_ + std::to_string(rank_id) + ".txt");
        this->WritePointData(point_outfile);
    }

    std::cout << " ---- Finish dividing " << this->CaseName() << " data ----"
              << "\n";
    return 0;
}

void DataPartitioner::InputMeshData(std::ifstream &infile) {
    for (int i = 0; i < 3; ++i) infile >> idimc[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; ++i) infile >> xyminw[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; ++i) infile >> xymaxw[i];
    infile.ignore(1000, '\n');

    infile >> nelem;
    for (int i = 0; i < 3; ++i) infile >> xyelemw[i];
    infile.ignore(1000, '\n');

    infile >> node;
    for (int i = 0; i < 3; ++i) infile >> xynodew[i];
    infile.ignore(1000, '\n');

    infile >> nodec;
    for (int i = 0; i < 3; ++i) infile >> xynodecw[i];
    infile.ignore(1000, '\n');
}

void DataPartitioner::OutputMeshData(std::ofstream &outfile, int rank_id) const {
    outfile << std::setw(15) << rank_id << "     --- myrank ---" << "\n";

    outfile << std::setw(15) << nelem << std::setw(15) << xyelem[0] << std::setw(15) << xyelem[1] << std::setw(15)
            << xyelem[2] << "     --- rank element ---" << "\n";

    outfile << std::setw(15) << node << std::setw(15) << xynode[0] << std::setw(15) << xynode[1] << std::setw(15)
            << xynode[2] << "     --- rank node ---" << "\n";

    outfile << std::setw(15) << nodec << std::setw(15) << xynodec[0] << std::setw(15) << xynodec[1] << std::setw(15)
            << xynodec[2] << "     --- rank control point ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << nxyr[i];
    outfile << "     --- rank number ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << xymin[i];
    outfile << "     --- rank minimum coordinate ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << xymax[i];
    outfile << "     --- rank maximum coordinate ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << aelemmin[i];
    outfile << "     --- rank minimum element ID ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << aelemmax[i];
    outfile << "     --- rank maximum element ID ---" << "\n";

    outfile << std::setw(15) << isb << "     --- overlap rank number ---" << "\n";
    if (isb != 0) {
        OutputVector(outfile, isb, naid);

        outfile << std::setw(15) << isubc << "     --- overlap control point number ---" << "\n";
        OutputVector(outfile, isb, nsbc);
        OutputVector(outfile, isubc, nsubc);
        OutputVector(outfile, nodec, dbc);

        outfile << std::setw(15) << isubl << "     --- overlap node number ---"
                << "\n";
        OutputVector(outfile, isb, nsbl);
        OutputVector(outfile, isubl, nsubl);
        OutputVector(outfile, node, dbl);
    }
}

void DataPartitioner::PartitionInitialDataset() {
    for (int i = 0; i < 3; ++i) {
        this->nexyr1_[i] = xyelemw[i] / nxyr[i];
        this->nexyr2_[i] = this->nexyr1_[i] + 1;
        this->nr2_[i] = xyelemw[i] % nxyr[i];
        this->nr1_[i] = nxyr[i] - this->nr2_[i];
        dxy[i] = (xymaxw[i] - xyminw[i]) / double(xyelemw[i]);
        this->drxy1_[i] = (xymaxw[i] - xyminw[i] - dxy[i] * double(this->nr2_[i])) / double(nxyr[i]);
        this->drxy2_[i] = this->drxy1_[i] + dxy[i];
        this->nxyr1_[i] = this->nexyr1_[i] + 1;
        this->nxyr2_[i] = this->nexyr2_[i] + 1;
        this->bound12_[i] = xyminw[i] + this->drxy1_[i] * double(this->nr1_[i]);
    }
}

void DataPartitioner::BcRenumber(BoundaryCondition &partition_bc, const BoundaryCondition &source_bc,
                                 const std::vector<int> &local_counts, const std::vector<int> &global_counts,
                                 const std::vector<int> &local_min, const std::vector<int> &local_max,
                                 bool compute_values) const {
    partition_bc.nbc.resize(source_bc.ibc);
    if (compute_values) { partition_bc.fbc.resize(source_bc.ibc); }

    std::vector<int> inxy(3), local_index(3);
    int nx = global_counts[0];
    int ny = global_counts[1];

    partition_bc.ibc = 0;
    for (int i = 0; i < source_bc.ibc; ++i) {
        int n = source_bc.nbc[i];
        inxy[2] = n / (nx * ny);
        inxy[1] = (n - inxy[2] * (nx * ny)) / nx;
        inxy[0] = n - inxy[2] * (nx * ny) - inxy[1] * nx;

        if (inxy[0] >= local_min[0] && inxy[0] <= local_max[0] && inxy[1] >= local_min[1] && inxy[1] <= local_max[1] &&
            inxy[2] >= local_min[2] && inxy[2] <= local_max[2]) {
            for (int j = 0; j < 3; ++j) { local_index[j] = inxy[j] - local_min[j]; }
            int local_id =
                local_counts[0] * local_counts[1] * local_index[2] + local_counts[0] * local_index[1] + local_index[0];
            partition_bc.nbc[partition_bc.ibc] = local_id;
            if (compute_values) { partition_bc.fbc[partition_bc.ibc] = source_bc.fbc[i]; }
            ++partition_bc.ibc;
        }
    }
}

void DataPartitioner::PointRenumber(int &local_num, std::vector<int> &global_id, //
                                    MaterialPoint &point, int rank_id) const {
    std::vector<int> irxy(3);
    global_id.resize(point.num);

    const int nrx = nxyr[0];
    const int nry = nxyr[1];

    local_num = 0;
    for (int ip = 0; ip < point.num; ++ip) {
        for (int i = 0; i < 3; ++i) {
            if (point.coord[ip][i] < this->bound12_[i]) {
                irxy[i] = floor((point.coord[ip][i] - xyminw[i]) / this->drxy1_[i]);
            } else {
                irxy[i] = this->nr1_[i] + floor((point.coord[ip][i] - this->bound12_[i]) / this->drxy2_[i]);
            }
        }
        int point_rank = nrx * nry * irxy[2] + nrx * irxy[1] + irxy[0];
        if (point_rank == rank_id) {
            global_id[local_num] = ip;
            ++local_num;
        }
    }
}

void DataPartitioner::MeshPartition(int rank_id) {
    const int nrx = nxyr[0];
    const int nry = nxyr[1];

    std::vector<int> irxy(3);
    irxy[2] = rank_id / (nrx * nry);
    irxy[1] = (rank_id - irxy[2] * (nrx * nry)) / nrx;
    irxy[0] = rank_id - irxy[2] * (nrx * nry) - irxy[1] * nrx;

    isb = 1;

    std::vector<std::vector<int>> ijkl(6, std::vector<int>(3));
    std::vector<std::vector<int>> ijkc(6, std::vector<int>(3));
    for (int ix = 0; ix < 3; ++ix) {
        if (irxy[ix] < this->nr1_[ix]) {
            xymin[ix] = xyminw[ix] + this->drxy1_[ix] * double(irxy[ix]);
            xymax[ix] = xymin[ix] + this->drxy1_[ix];
            aelemmin[ix] = this->nexyr1_[ix] * irxy[ix];
            aelemmax[ix] = aelemmin[ix] + this->nexyr1_[ix] - 1;
            xyelem[ix] = this->nexyr1_[ix];
        } else {
            xymin[ix] = this->bound12_[ix] + this->drxy2_[ix] * double(irxy[ix] - this->nr1_[ix]);
            xymax[ix] = xymin[ix] + this->drxy2_[ix];
            aelemmin[ix] = this->nexyr1_[ix] * this->nr1_[ix] + this->nexyr2_[ix] * (irxy[ix] - this->nr1_[ix]);
            aelemmax[ix] = aelemmin[ix] + this->nexyr2_[ix] - 1;
            xyelem[ix] = this->nexyr2_[ix];
        }

        xynode[ix] = xyelem[ix] + 1;
        xynodec[ix] = xyelem[ix] + idimc[ix];

        if (nxyr[ix] == 1) {
            int temp_ijkl[] = {-1, -1, 0, xynode[ix] - 1, -1, -1};
            int temp_ijkc[] = {-1, -1, 0, xynodec[ix] - 1, -1, -1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
        } else if (irxy[ix] == 0) {
            int temp_ijkl[] = {-1, -1, 0, xynode[ix] - 1, xynode[ix] - 1, xynode[ix] - 1};
            int temp_ijkc[] = {-1, -1, 0, xynodec[ix] - 1, xynodec[ix] - idimc[ix], xynodec[ix] - 1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            isb *= 2;
        } else if (irxy[ix] == nxyr[ix] - 1) {
            int temp_ijkl[] = {0, 0, 0, xynode[ix] - 1, -1, -1};
            int temp_ijkc[] = {0, idimc[ix] - 1, 0, xynodec[ix] - 1, -1, -1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            isb *= 2;
        } else {
            int temp_ijkl[] = {0, 0, 0, xynode[ix] - 1, xynode[ix] - 1, xynode[ix] - 1};
            int temp_ijkc[] = {0, idimc[ix] - 1, 0, xynodec[ix] - 1, xynodec[ix] - idimc[ix], xynodec[ix] - 1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            isb *= 3;
        }

        this->inxyminl_[ix] = aelemmin[ix];
        this->inxymaxl_[ix] = aelemmax[ix] + 1;
        this->inxyminc_[ix] = aelemmin[ix];
        this->inxymaxc_[ix] = aelemmax[ix] + idimc[ix];
    }

    node = xynode[0] * xynode[1] * xynode[2];
    nodec = xynodec[0] * xynodec[1] * xynodec[2];
    nelem = xyelem[0] * xyelem[1] * xyelem[2];

    isb--;
    naid.resize(isb);
    nsbc.resize(isb);
    nsbl.resize(isb);
    nsubc.resize(nodec * isb);
    nsubl.resize(node * isb);
    VectorAssign(nodec, dbc, 1.0e0);
    VectorAssign(node, dbl, 1.0e0);

    isubc = 0;
    isubl = 0;
    isb = 0;
    for (int kk = 0; kk < 3; ++kk) {
        int ksc = ijkc[2 * kk][2];
        int kec = ijkc[2 * kk + 1][2];
        int ksl = ijkl[2 * kk][2];
        int kel = ijkl[2 * kk + 1][2];
        if (ksc == -1) continue;
        for (int jj = 0; jj < 3; ++jj) {
            int jsc = ijkc[2 * jj][1];
            int jec = ijkc[2 * jj + 1][1];
            int jsl = ijkl[2 * jj][1];
            int jel = ijkl[2 * jj + 1][1];
            if (jsc == -1) continue;
            for (int ii = 0; ii < 3; ++ii) {
                if (ii * jj * kk == 1) continue;
                int isc = ijkc[2 * ii][0];
                int iec = ijkc[2 * ii + 1][0];
                int isl = ijkl[2 * ii][0];
                int iel = ijkl[2 * ii + 1][0];
                if (isc == -1) continue;

                naid[isb] = rank_id + (ii - 1) + (jj - 1) * nxyr[0] + (kk - 1) * nxyr[0] * nxyr[1];
                nsbc[isb] = (iec - isc + 1) * (jec - jsc + 1) * (kec - ksc + 1);
                nsbl[isb] = (iel - isl + 1) * (jel - jsl + 1) * (kel - ksl + 1);

                for (int k = ksl; k <= kel; ++k) {
                    for (int j = jsl; j <= jel; ++j) {
                        for (int i = isl; i <= iel; ++i) {
                            int nn = i + j * xynode[0] + k * xynode[0] * xynode[1];
                            nsubl[isubl++] = nn;
                            dbl[nn] += 1.0e0;
                        }
                    }
                }

                for (int k = ksc; k <= kec; ++k) {
                    for (int j = jsc; j <= jec; ++j) {
                        for (int i = isc; i <= iec; ++i) {
                            int nn = i + j * xynodec[0] + k * xynodec[0] * xynodec[1];
                            nsubc[isubc++] = nn;
                            dbc[nn] += 1.0e0;
                        }
                    }
                }
                ++isb;
            }
        }
    }
}

void DataPartitioner::LoadPartitionParameters(std::ifstream &infile) {
    for (int i = 0; i < 3; ++i) { infile >> nxyr[i]; }
    infile >> this->grid_input_file_;
    infile >> this->point_input_file_;
    infile >> this->grid_output_prefix_;
    infile >> this->point_output_prefix_;
}

void DataPartitioner::CreateOutputDirectories() const {
    const std::filesystem::path grid_parent = std::filesystem::path(this->grid_output_prefix_).parent_path();
    if (!grid_parent.empty()) { std::filesystem::create_directories(grid_parent); }

    const std::filesystem::path point_parent = std::filesystem::path(this->point_output_prefix_).parent_path();
    if (!point_parent.empty()) { std::filesystem::create_directories(point_parent); }
}
