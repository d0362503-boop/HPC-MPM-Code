#include "data/divide/data_partitioner.h"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

int DataPartitioner::Run() {
    Partition &part = this->grid_.part;

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

    const int rank_count = part.nxyr[0] * part.nxyr[1] * part.nxyr[2];
    for (int rank_id = 0; rank_id < rank_count; ++rank_id) {
        this->PartitionProcess(rank_id);

        std::ofstream grid_outfile = OpenOutputFile(this->grid_output_prefix_ + std::to_string(rank_id) + ".txt");
        this->OutputMeshData(grid_outfile, rank_id);
        this->WriteBoundaryData(grid_outfile);

        std::ofstream point_outfile = OpenOutputFile(this->point_output_prefix_ + std::to_string(rank_id) + ".txt");
        this->WritePointData(point_outfile);
    }

    std::cout << " ---- Finish dividing " << this->CaseName() << " data ----"
              << "\n";

    return 0;
}

void DataPartitioner::InputMeshData(std::ifstream &infile) {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;

    for (int i = 0; i < 3; ++i) infile >> cmesh.idimc[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; ++i) infile >> mesh.xyminw[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; ++i) infile >> mesh.xymaxw[i];
    infile.ignore(1000, '\n');

    infile >> mesh.nelem;
    for (int i = 0; i < 3; ++i) infile >> mesh.xyelemw[i];
    infile.ignore(1000, '\n');

    infile >> mesh.node;
    for (int i = 0; i < 3; ++i) infile >> mesh.xynodew[i];
    infile.ignore(1000, '\n');

    infile >> cmesh.nodec;
    for (int i = 0; i < 3; ++i) infile >> cmesh.xynodecw[i];
    infile.ignore(1000, '\n');

    return;
}

void DataPartitioner::OutputMeshData(std::ofstream &outfile, int rank_id) const {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;
    Partition &part = this->grid_.part;

    outfile << std::setw(15) << rank_id << "     --- myrank ---" << "\n";

    outfile << std::setw(15) << mesh.nelem << std::setw(15) << mesh.xyelem[0] << std::setw(15) << mesh.xyelem[1] << std::setw(15)
            << mesh.xyelem[2] << "     --- rank element ---" << "\n";

    outfile << std::setw(15) << mesh.node << std::setw(15) << mesh.xynode[0] << std::setw(15) << mesh.xynode[1] << std::setw(15)
            << mesh.xynode[2] << "     --- rank node ---" << "\n";

    outfile << std::setw(15) << cmesh.nodec << std::setw(15) << cmesh.xynodec[0] << std::setw(15) << cmesh.xynodec[1] << std::setw(15)
            << cmesh.xynodec[2] << "     --- rank control point ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << part.nxyr[i];
    outfile << "     --- rank number ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << mesh.xymin[i];
    outfile << "     --- rank minimum coordinate ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << mesh.xymax[i];
    outfile << "     --- rank maximum coordinate ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << mesh.aelemmin[i];
    outfile << "     --- rank minimum element ID ---" << "\n";

    for (int i = 0; i < 3; ++i) outfile << std::setw(15) << mesh.aelemmax[i];
    outfile << "     --- rank maximum element ID ---" << "\n";

    outfile << std::setw(15) << part.isb << "     --- overlap rank number ---" << "\n";
    if (part.isb != 0) {
        OutputVector(outfile, part.isb, part.naid);

        outfile << std::setw(15) << part.isubc << "     --- overlap control point number ---" << "\n";
        OutputVector(outfile, part.isb, part.nsbc);
        OutputVector(outfile, part.isubc, part.nsubc);
        OutputVector(outfile, cmesh.nodec, part.dbc);

        outfile << std::setw(15) << part.isubl << "     --- overlap node number ---"
                << "\n";
        OutputVector(outfile, part.isb, part.nsbl);
        OutputVector(outfile, part.isubl, part.nsubl);
        OutputVector(outfile, mesh.node, part.dbl);
    }

    return;
}

void DataPartitioner::PartitionInitialDataset() {
    Mesh &mesh = this->grid_.mesh;
    Partition &part = this->grid_.part;

    for (int i = 0; i < 3; ++i) {
        this->nexyr1_[i] = mesh.xyelemw[i] / part.nxyr[i];
        this->nexyr2_[i] = this->nexyr1_[i] + 1;
        this->nr2_[i] = mesh.xyelemw[i] % part.nxyr[i];
        this->nr1_[i] = part.nxyr[i] - this->nr2_[i];
        mesh.dxy[i] = (mesh.xymaxw[i] - mesh.xyminw[i]) / double(mesh.xyelemw[i]);
        this->drxy1_[i] = (mesh.xymaxw[i] - mesh.xyminw[i] - mesh.dxy[i] * double(this->nr2_[i])) / double(part.nxyr[i]);
        this->drxy2_[i] = this->drxy1_[i] + mesh.dxy[i];
        this->nxyr1_[i] = this->nexyr1_[i] + 1;
        this->nxyr2_[i] = this->nexyr2_[i] + 1;
        this->bound12_[i] = mesh.xyminw[i] + this->drxy1_[i] * double(this->nr1_[i]);
    }

    return;
}

void DataPartitioner::BCRenumber(BoundaryCondition &partition_bc, const BoundaryCondition &source_bc,
                                 const std::array<int, 3> &local_counts, const std::array<int, 3> &global_counts,
                                 const std::array<int, 3> &local_min, const std::array<int, 3> &local_max,
                                 bool compute_values) {

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

        if (inxy[0] >= local_min[0] && inxy[0] <= local_max[0] && //
            inxy[1] >= local_min[1] && inxy[1] <= local_max[1] && //
            inxy[2] >= local_min[2] && inxy[2] <= local_max[2]) {
            for (int j = 0; j < 3; ++j) { local_index[j] = inxy[j] - local_min[j]; }
            int local_id =
                local_counts[0] * local_counts[1] * local_index[2] + local_counts[0] * local_index[1] + local_index[0];
            partition_bc.nbc[partition_bc.ibc] = local_id;
            if (compute_values) { partition_bc.fbc[partition_bc.ibc] = source_bc.fbc[i]; }
            ++partition_bc.ibc;
        }
    }

    return;
}

void DataPartitioner::PointRenumber(int &local_num, std::vector<int> &global_id, //
                                    MaterialPoint &point, int rank_id) const {
    Mesh &mesh = this->grid_.mesh;
    Partition &part = this->grid_.part;

    std::vector<int> irxy(3);
    global_id.resize(point.num);

    const int nrx = part.nxyr[0];
    const int nry = part.nxyr[1];

    local_num = 0;
    for (int ip = 0; ip < point.num; ++ip) {
        for (int i = 0; i < 3; ++i) {
            if (point.coord[ip][i] < this->bound12_[i]) {
                irxy[i] = floor((point.coord[ip][i] - mesh.xyminw[i]) / this->drxy1_[i]);
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

    return;
}

void DataPartitioner::MeshPartition(int rank_id) {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;
    Partition &part = this->grid_.part;

    const int nrx = part.nxyr[0];
    const int nry = part.nxyr[1];

    std::vector<int> irxy(3);
    irxy[2] = rank_id / (nrx * nry);
    irxy[1] = (rank_id - irxy[2] * (nrx * nry)) / nrx;
    irxy[0] = rank_id - irxy[2] * (nrx * nry) - irxy[1] * nrx;

    part.isb = 1;

    std::vector<std::vector<int>> ijkl(6, std::vector<int>(3));
    std::vector<std::vector<int>> ijkc(6, std::vector<int>(3));
    for (int ix = 0; ix < 3; ++ix) {
        if (irxy[ix] < this->nr1_[ix]) {
            mesh.xymin[ix] = mesh.xyminw[ix] + this->drxy1_[ix] * double(irxy[ix]);
            mesh.xymax[ix] = mesh.xymin[ix] + this->drxy1_[ix];
            mesh.aelemmin[ix] = this->nexyr1_[ix] * irxy[ix];
            mesh.aelemmax[ix] = mesh.aelemmin[ix] + this->nexyr1_[ix] - 1;
            mesh.xyelem[ix] = this->nexyr1_[ix];
        } else {
            mesh.xymin[ix] = this->bound12_[ix] + this->drxy2_[ix] * double(irxy[ix] - this->nr1_[ix]);
            mesh.xymax[ix] = mesh.xymin[ix] + this->drxy2_[ix];
            mesh.aelemmin[ix] = this->nexyr1_[ix] * this->nr1_[ix] + this->nexyr2_[ix] * (irxy[ix] - this->nr1_[ix]);
            mesh.aelemmax[ix] = mesh.aelemmin[ix] + this->nexyr2_[ix] - 1;
            mesh.xyelem[ix] = this->nexyr2_[ix];
        }

        mesh.xynode[ix] = mesh.xyelem[ix] + 1;
        cmesh.xynodec[ix] = mesh.xyelem[ix] + cmesh.idimc[ix];

        if (part.nxyr[ix] == 1) {
            int temp_ijkl[] = {-1, -1, 0, mesh.xynode[ix] - 1, -1, -1};
            int temp_ijkc[] = {-1, -1, 0, cmesh.xynodec[ix] - 1, -1, -1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
        } else if (irxy[ix] == 0) {
            int temp_ijkl[] = {-1, -1, 0, mesh.xynode[ix] - 1, mesh.xynode[ix] - 1, mesh.xynode[ix] - 1};
            int temp_ijkc[] = {-1, -1, 0, cmesh.xynodec[ix] - 1, cmesh.xynodec[ix] - cmesh.idimc[ix], cmesh.xynodec[ix] - 1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            part.isb *= 2;
        } else if (irxy[ix] == part.nxyr[ix] - 1) {
            int temp_ijkl[] = {0, 0, 0, mesh.xynode[ix] - 1, -1, -1};
            int temp_ijkc[] = {0, cmesh.idimc[ix] - 1, 0, cmesh.xynodec[ix] - 1, -1, -1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            part.isb *= 2;
        } else {
            int temp_ijkl[] = {0, 0, 0, mesh.xynode[ix] - 1, mesh.xynode[ix] - 1, mesh.xynode[ix] - 1};
            int temp_ijkc[] = {0, cmesh.idimc[ix] - 1, 0, cmesh.xynodec[ix] - 1, cmesh.xynodec[ix] - cmesh.idimc[ix], cmesh.xynodec[ix] - 1};
            for (int i = 0; i < 6; ++i) {
                ijkc[i][ix] = temp_ijkc[i];
                ijkl[i][ix] = temp_ijkl[i];
            }
            part.isb *= 3;
        }

        this->inxyminl_[ix] = mesh.aelemmin[ix];
        this->inxymaxl_[ix] = mesh.aelemmax[ix] + 1;
        this->inxyminc_[ix] = mesh.aelemmin[ix];
        this->inxymaxc_[ix] = mesh.aelemmax[ix] + cmesh.idimc[ix];
    }

    mesh.node = mesh.xynode[0] * mesh.xynode[1] * mesh.xynode[2];
    cmesh.nodec = cmesh.xynodec[0] * cmesh.xynodec[1] * cmesh.xynodec[2];
    mesh.nelem = mesh.xyelem[0] * mesh.xyelem[1] * mesh.xyelem[2];

    part.isb--;
    part.naid.resize(part.isb);
    part.nsbc.resize(part.isb);
    part.nsbl.resize(part.isb);
    part.nsubc.resize(cmesh.nodec * part.isb);
    part.nsubl.resize(mesh.node * part.isb);
    VectorAssign(cmesh.nodec, part.dbc, 1.0e0);
    VectorAssign(mesh.node, part.dbl, 1.0e0);

    part.isubc = 0;
    part.isubl = 0;
    part.isb = 0;
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

                part.naid[part.isb] = rank_id + (ii - 1) + (jj - 1) * part.nxyr[0] + (kk - 1) * part.nxyr[0] * part.nxyr[1];
                part.nsbc[part.isb] = (iec - isc + 1) * (jec - jsc + 1) * (kec - ksc + 1);
                part.nsbl[part.isb] = (iel - isl + 1) * (jel - jsl + 1) * (kel - ksl + 1);

                for (int k = ksl; k <= kel; ++k) {
                    for (int j = jsl; j <= jel; ++j) {
                        for (int i = isl; i <= iel; ++i) {
                            int nn = i + j * mesh.xynode[0] + k * mesh.xynode[0] * mesh.xynode[1];
                            part.nsubl[part.isubl++] = nn;
                            part.dbl[nn] += 1.0e0;
                        }
                    }
                }

                for (int k = ksc; k <= kec; ++k) {
                    for (int j = jsc; j <= jec; ++j) {
                        for (int i = isc; i <= iec; ++i) {
                            int nn = i + j * cmesh.xynodec[0] + k * cmesh.xynodec[0] * cmesh.xynodec[1];
                            part.nsubc[part.isubc++] = nn;
                            part.dbc[nn] += 1.0e0;
                        }
                    }
                }
                ++part.isb;
            }
        }
    }

    return;
}

void DataPartitioner::LoadPartitionParameters(std::ifstream &infile) {
    Partition &part = this->grid_.part;

    for (int i = 0; i < 3; ++i) { infile >> part.nxyr[i]; }
    infile >> this->grid_input_file_;
    infile >> this->point_input_file_;
    infile >> this->grid_output_prefix_;
    infile >> this->point_output_prefix_;

    return;
}

void DataPartitioner::CreateOutputDirectories() const {

    const std::filesystem::path grid_parent = std::filesystem::path(this->grid_output_prefix_).parent_path();
    if (!grid_parent.empty()) { std::filesystem::create_directories(grid_parent); }

    const std::filesystem::path point_parent = std::filesystem::path(this->point_output_prefix_).parent_path();
    if (!point_parent.empty()) { std::filesystem::create_directories(point_parent); }

    return;
}
