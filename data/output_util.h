#ifndef DATA_OUTPUT_UTIL_H_
#define DATA_OUTPUT_UTIL_H_

#include <fstream>
#include <string>

#include "../module/material_point.h"

void WriteGlobalMeshHeader(std::ofstream &outfile);

#ifdef HAVE_HDF5
void WriteVtkHdf5Mesh(const std::string &filename);

void WriteVtkHdf5Points(const std::string &filename, const MaterialPoint &point);
#endif

#endif // DATA_OUTPUT_UTIL_H_
