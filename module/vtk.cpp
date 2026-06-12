#include <iostream>
#include <fstream>
#include <vector>
#include <arpa/inet.h> 
#include <iomanip>
#include "vtk.h"
#include "mesh.h"
#include "material_point.h"

using namespace std;

void OutputVtkMeshInitialize(ofstream& vtkfile) {

    vtkfile << "# vtk DataFile Version 3.0" << "\n";
    vtkfile << "Unstructured Grid for Seminor" << "\n";
    vtkfile << "BINARY" << "\n";

    vtkfile << "DATASET UNSTRUCTURED_GRID" << "\n";
    vtkfile << "POINTS " << node << " float" << "\n";
    for (int j = 0; j < node; ++j) {
        float coords[3] = {static_cast<float>(xyn[j][0]), static_cast<float>(xyn[j][1]), 
                           static_cast<float>(xyn[j][2])};
        for (int i = 0; i < 3; ++i) {
            uint32_t val = htonl(*reinterpret_cast<uint32_t*>(&coords[i]));  // For VTK file BIG_ENDIAN
            vtkfile.write(reinterpret_cast<char*>(&val), sizeof(val)); 
        }
    }
    vtkfile << "\n";

    vtkfile << "CELLS " << nelem << " " << nelem * 9 << "\n";
    for (int j = 0; j < nelem; ++j) {
        uint32_t num = htonl(8);  // Eight nodes in each element
        vtkfile.write(reinterpret_cast<char*>(&num), sizeof(num));
        for (int i = 0; i < 8; ++i) {
            uint32_t idx = htonl(nc[j][i]);  // VTK ID start form 0
            vtkfile.write(reinterpret_cast<char*>(&idx), sizeof(idx));
        }
    }
    vtkfile << "\n";

    vtkfile << "CELL_TYPES " << nelem << "\n";
    for (int j = 0; j < nelem; ++j) {
        uint32_t type = htonl(12);  // VTK_HEXAHEDRON
        vtkfile.write(reinterpret_cast<char*>(&type), sizeof(type));
    }
    vtkfile << "\n";

    return;
}

void OutputVtkPointInitialize(ofstream& vtkfile, MaterialPoint& point) {

    vtkfile << "# vtk DataFile Version 3.0" << "\n";
    vtkfile << "Unstructured Grid for Seminor" << "\n";
    vtkfile << "BINARY" << "\n";

    vtkfile << "DATASET UNSTRUCTURED_GRID" << "\n";
    vtkfile << "POINTS " << point.num << " float\n";
    OutputVtkVector(vtkfile, point.num, point.coord);

    vtkfile << "CELLS " << point.num << " " << point.num * 2 << "\n";
    for (int j = 0; j < point.num; ++j) {
        uint32_t num = htonl(1);
        vtkfile.write(reinterpret_cast<char*>(&num), sizeof(num));
        uint32_t idx = htonl(j);
        vtkfile.write(reinterpret_cast<char*>(&idx), sizeof(idx));
    }
    vtkfile << "\n";

    vtkfile << "CELL_TYPES " << point.num << "\n";
    for (int j = 0; j < point.num; ++j) {
        uint32_t type = htonl(1); // VTK_VERTEX
        vtkfile.write(reinterpret_cast<char*>(&type), sizeof(type));
    }
    vtkfile << "\n";

    return;
}
