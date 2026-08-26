#include "data/generate/data_generator.h"

#include <iostream>
#include <stdexcept>

#include "data/generate/output_util.h"
#include "module/data_io.h"
#include "module/dataset.h"
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

    this->InitializeGridGeometry();

    BuildMesh();

    BuildControlPoint();

    this->CreateBCs();

    this->CreateParticles();

    return;
}

void DataGenerator::CreateBCs() { throw std::logic_error("DataGenerator requires a CreateBCs override."); }

void DataGenerator::CreateParticles() { throw std::logic_error("DataGenerator requires a CreateParticles override."); }

void DataGenerator::InitializeGridGeometry() {

    for (int i = 0; i < 3; ++i) {
        dxy[i] = (xymax[i] - xymin[i]) / static_cast<double>(xyelem[i]);
        aelemmin[i] = 0;
        aelemmax[i] = xyelem[i] - 1;
        xynode[i] = xyelem[i] + 1;
        xynodec[i] = xyelem[i] + idimc[i];
    }
    nelem = xyelem[0] * xyelem[1] * xyelem[2];
    node = xynode[0] * xynode[1] * xynode[2];
    nodec = xynodec[0] * xynodec[1] * xynodec[2];

    return;
}

void DataGenerator::WriteGridDataFile(const std::string &filename) {

    std::ofstream outfile = OpenOutputFile(filename);
    WriteGlobalMeshHeader(outfile);
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
