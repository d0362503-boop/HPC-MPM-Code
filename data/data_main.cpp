#include "../module/bc.h"
#include "../module/dataset.h"
#include "../module/material_point.h"
#include "../module/mesh.h"
#include "../module/solid/implicit/implicit_mpm_solid.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace implicitmpm;

void create_data();
void gloabl_mesh_data_output(string filename);
void mesh_vtk(string filename);
void point_vtk(string filename, MaterialPoint &point);

int main() {

    cout << " ---- Start making data ----" << "\n";

    create_data();

    {
        cout << "Making grid file" << "\n";
        // string gridfile = "fluid/griddata.txt";     // Grid output filename
        // string gridfile = "solid/griddata.txt";     // Grid output filename
        string gridfile = "fsi/griddata.txt"; // Grid output filename
        gloabl_mesh_data_output(gridfile);
    }

    {
        cout << "Making VTK file" << "\n";
        // string gridfile = "fluid/grid.vtk";     // Grid output VTK filename
        // string gridfile = "solid/grid.vtk";     // Grid output VTK filename
        string gridfile = "fsi/grid.vtk"; // Grid output VTK filename
        mesh_vtk(gridfile);

        string pointfile;
        // pointfile = "solid/sp.vtk";      // Solid point output VTK filename
        pointfile = "fsi/sp.vtk"; // Solid point output VTK filename
        point_vtk(pointfile, sp);

        // pointfile = "fluid/wp.vtk";      // Water point output VTK filename
        // point_vtk (pointfile, wp);
    }

    cout << " ---- Finish making data ----" << "\n";

    return 0;
}

void global_bc_data_output(ofstream &outfile);

void gloabl_mesh_data_output(string filename) {

    ofstream gridfile;
    gridfile.open(filename);

    gridfile.flags(ios::right | ios::scientific); // --- output format ---
    for (int i = 0; i < 3; ++i) gridfile << setw(10) << idimc[i];
    gridfile << "\n";

    for (int i = 0; i < 3; ++i) gridfile << setw(15) << xymin[i];
    gridfile << "\n";

    for (int i = 0; i < 3; ++i) gridfile << setw(15) << xymax[i];
    gridfile << "\n";

    gridfile << setw(10) << nelem;
    for (int i = 0; i < 3; ++i) gridfile << setw(10) << xyelem[i];
    gridfile << "\n";

    gridfile << setw(10) << node;
    for (int i = 0; i < 3; ++i) gridfile << setw(10) << xynode[i];
    gridfile << "\n";

    gridfile << setw(10) << nodec;
    for (int i = 0; i < 3; ++i) gridfile << setw(10) << xynodec[i];
    gridfile << "\n";

    global_bc_data_output(gridfile);

    gridfile.close();

    return;
}

void mesh_vtk(string filename) {

    ofstream vtkfile;
    vtkfile.open(filename);

    OutputVtkMeshInitialize(vtkfile);

    vtkfile.close();

    return;
}

void point_vtk(string filename, MaterialPoint &point) {

    ofstream vtkfile;
    vtkfile.open(filename);

    OutputVtkPointInitialize(vtkfile, point);

    vtkfile << "POINT_DATA " << point.num << "\n";
    vtkfile << "SCALARS matid float" << "\n";
    OutputVtkScalar(vtkfile, point.num, point.matid);

    if (point.surf_point.empty() == false) {
        vtkfile << "SCALARS surf_flag float" << "\n";
        OutputVtkScalar(vtkfile, point.num, point.surf_point);
    }

    vtkfile.close();

    return;
}
