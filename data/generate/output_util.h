#ifndef DATA_OUTPUT_UTIL_H_
#define DATA_OUTPUT_UTIL_H_

#include <fstream>
#include <string>

#include "module/material_point.h"

/**
 * @brief Write the global mesh header (B-spline orders, box, element/node/control-point counts).
 * @param outfile Output stream of the generated grid-data file.
 * @param grid    Generated global grid.
 */
void WriteGlobalMeshHeader(std::ofstream &outfile, const Grid &grid);

#ifdef HAVE_HDF5
/**
 * @brief Write the background-mesh topology to a VTK HDF5 file.
 * @param filename Output file name.
 * @param mesh     Generated background mesh (node coordinates and connectivity).
 */
void WriteVTKHDFMesh(const std::string &filename, const Mesh &mesh);

void WriteVTKHDFPoints(const std::string &filename, const MaterialPoint &point);
#endif

#endif // DATA_OUTPUT_UTIL_H_
