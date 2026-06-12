#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <arpa/inet.h> 
#include <iomanip>
#include <string>
#include "material_point.h"

using namespace std;

void OutputVtkMeshInitialize(ofstream& vtkfile);

void OutputVtkPointInitialize(ofstream& vtkfile, MaterialPoint& point);

template <typename T>
void OutputVtkScalar (ofstream& vtkfile, int num, vector<T>& variable) {

    vtkfile << "LOOKUP_TABLE default \n";
    for (int i = 0; i < num; ++i) {
        float val_float = static_cast<float>(variable[i]);
        uint32_t val = htonl(*reinterpret_cast<uint32_t*>(&val_float));
        vtkfile.write(reinterpret_cast<char*>(&val), sizeof(val));
    }
    vtkfile << "\n";

    return;
}

template <typename T>
void OutputVtkVector (ofstream& vtkfile, int num, vector<array<T, 3>>& variable) {

    for (int j = 0; j < num; ++j) {
        float vec[3] = {static_cast<float>(variable[j][0]), static_cast<float>(variable[j][1]),
                        static_cast<float>(variable[j][2])};
        for (int i = 0; i < 3; ++i) {
            uint32_t val = htonl(*reinterpret_cast<uint32_t*>(&vec[i]));
            vtkfile.write(reinterpret_cast<char*>(&val), sizeof(val));
        }
    }
    vtkfile << "\n";

    return;
}


