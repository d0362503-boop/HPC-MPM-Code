#include "data/generate/data_generator.h"

#include <iostream>
#include <stdexcept>

#include "data/generate/output_util.h"
#include "module/data_io.h"
#include "module/mesh.h"

int DataGenerator::Run() {

    std::cout << " ---- Start making " << this->CaseName() << " data ----"
              << "\n";

    this->LoadInput();

    this->BuildData();

    this->WriteTextOutputs();

    this->WriteVisualizationOutputs();

    std::cout << " ---- Finish making " << this->CaseName() << " data ----"
              << "\n";

    return 0;
}

void DataGenerator::BuildData() {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;

    this->InitializeGridGeometry();

    BuildMesh(mesh);

    BuildControlPoint(mesh, cmesh);

    this->CreateBCs();

    this->CreateParticles();

    return;
}

void DataGenerator::CreateBCs() { throw std::logic_error("DataGenerator requires a CreateBCs override."); }

void DataGenerator::CreateParticles() { throw std::logic_error("DataGenerator requires a CreateParticles override."); }

void DataGenerator::InitializeGridGeometry() {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;

    for (int i = 0; i < 3; ++i) {
        mesh.dxy[i] = (mesh.xymax[i] - mesh.xymin[i]) / static_cast<double>(mesh.xyelem[i]);
        mesh.aelemmin[i] = 0;
        mesh.aelemmax[i] = mesh.xyelem[i] - 1;
        mesh.xynode[i] = mesh.xyelem[i] + 1;
        cmesh.xynodec[i] = mesh.xyelem[i] + cmesh.idimc[i];
    }
    mesh.nelem = mesh.xyelem[0] * mesh.xyelem[1] * mesh.xyelem[2];
    mesh.node = mesh.xynode[0] * mesh.xynode[1] * mesh.xynode[2];
    cmesh.nodec = cmesh.xynodec[0] * cmesh.xynodec[1] * cmesh.xynodec[2];

    return;
}

void DataGenerator::WriteGridDataFile(const std::string &filename) {

    std::ofstream outfile = OpenOutputFile(filename);
    WriteGlobalMeshHeader(outfile, this->grid_);
    this->WriteBCData(outfile);

    return;
}

void DataGenerator::WriteTextOutputs() {

    std::cout << "Making grid file" << "\n";
    this->WriteGridDataFile("griddata.txt");

    std::cout << "Making point file" << "\n";
    std::ofstream pointfile = OpenOutputFile("pointdata.txt");
    this->WritePointData(pointfile);

    return;
}

void DataGenerator::WritePointData(std::ofstream &) {

    throw std::logic_error("DataGenerator requires a WritePointData override.");
}
