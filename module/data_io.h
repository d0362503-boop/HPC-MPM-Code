#pragma once

#include "material_point.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

/**
 * @brief Open an input file and verify it was successfully opened.
 * @param filename Path to the input file.
 * @return Opened input file stream.
 * @throws std::runtime_error if the file cannot be opened.
 */
std::ifstream OpenInputFile(const std::string &filename);

/**
 * @brief Open an output file and verify it was successfully opened.
 * @param filename Path to the output file.
 * @return Opened output file stream.
 * @throws std::runtime_error if the file cannot be opened.
 */
std::ofstream OpenOutputFile(const std::string &filename);

/**
 * @brief Read a vector of scalars from an input stream.
 * @param infile Input file stream.
 * @param size0  Number of entries to read.
 * @param variable Output vector (must be pre-sized to at least `size0`).
 */
template <typename T> void InputVector(std::ifstream &infile, const int size0, std::vector<T> &variable) {

    for (int i = 0; i < size0; i++) {
        infile >> variable[i];
        infile.ignore(1000, '\n');
    }
    return;
}

/**
 * @brief Read a vector of fixed-size arrays from an input stream.
 * @param infile   Input file stream.
 * @param size0    Number of array entries to read.
 * @param variable Output vector of arrays (must be pre-sized to at least `size0`).
 */
template <typename T, std::size_t N>
void InputVector(std::ifstream &infile, const int size0, std::vector<std::array<T, N>> &variable) {

    for (int i = 0; i < size0; i++) {
        for (int j = 0; j < N; j++) {
            if constexpr (std::is_arithmetic_v<T>) {
                infile >> variable[i][j];
            } else {
                for (size_t k = 0; k < variable[i][j].size(); k++) { infile >> variable[i][j][k]; }
            }
        }
        infile.ignore(1000, '\n');
    }
    return;
}

/**
 * @brief Write a vector of scalars to an output stream.
 * @param outfile Output file stream.
 * @param size0   Number of entries to write.
 * @param variable Vector of values to write.
 */
template <typename T> void OutputVector(std::ofstream &outfile, const int size0, std::vector<T> &variable) {

    for (int i = 0; i < size0; i++) { outfile << std::setw(15) << variable[i] << "\n"; }
    return;
}

/**
 * @brief Write a vector of fixed-size arrays to an output stream.
 * @param outfile Output file stream.
 * @param size0   Number of array entries to write.
 * @param variable Vector of arrays to write.
 */
template <typename T, std::size_t N>
void OutputVector(std::ofstream &outfile, const int size0, std::vector<std::array<T, N>> &variable) {

    for (int i = 0; i < size0; i++) {
        for (int j = 0; j < N; j++) {
            if constexpr (std::is_arithmetic_v<T>) {
                outfile << std::setw(15) << variable[i][j];
            } else {
                for (size_t k = 0; k < variable[i][j].size(); k++) { outfile << std::setw(15) << variable[i][j][k]; }
            }
        }
        outfile << "\n";
    }
    return;
}

/**
 * @brief Read per-rank MPI domain-decomposition data from a partitioned grid file.
 * @param infile Input file stream positioned at the grid-data header.
 */
void InputParaGriddata(std::ifstream &infile);

/**
 * @brief Print a progress message to stdout showing current step and view index.
 * @param iview Output view index.
 * @param istep Current time step.
 */
inline void OutputMessage(int iview, int istep) {

    std::cout << " ===== Now Computing ===== " << "\n";
    std::cout << "istep:" << std::setw(10) << istep << "\n";
    std::cout << "time:" << std::setw(15) << real_time << std::setw(10) //
              << "iview:" << std::setw(10) << iview << "\n";
    std::cout << "========================== " << "\n";

    return;
}
