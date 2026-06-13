#include "data_generator.h"

#include <iostream>
#include <stdexcept>

#include "../module/dataset.h"
#include "output_util.h"

int DataGenerator::Run() {
    std::cout << " ---- Start making " << this->CaseName() << " data ----"
              << "\n";

    this->Initialize();
    this->LoadInput();
    this->BuildData();
    this->WriteTextOutputs();
    this->WriteVisualizationOutputs();
    this->Finalize();

    std::cout << " ---- Finish making " << this->CaseName() << " data ----"
              << "\n";
    return 0;
}

std::ifstream DataGenerator::OpenInputFile(const std::string &filename) const {
    std::ifstream infile(filename);
    if (!infile.is_open()) { throw std::runtime_error("Failed to open input file: " + filename); }
    return infile;
}

std::ofstream DataGenerator::OpenOutputFile(const std::string &filename) const {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) { throw std::runtime_error("Failed to open output file: " + filename); }
    outfile.flags(std::ios::right | std::ios::scientific);
    return outfile;
}

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
}

void DataGenerator::WriteGridDataFile(const std::string &filename) {
    std::ofstream outfile = this->OpenOutputFile(filename);
    WriteGlobalMeshHeader(outfile);
    this->WriteBcData(outfile);
}
