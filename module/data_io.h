#pragma once

#include "material_point.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

std::ifstream OpenInputFile(const std::string &filename);

std::ofstream OpenOutputFile(const std::string &filename);

template <typename T> void InputVector(std::ifstream &infile, const int size0, std::vector<T> &variable) {

    for (int i = 0; i < size0; i++) {
        infile >> variable[i];
        infile.ignore(1000, '\n');
    }
    return;
}

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

template <typename T> void OutputVector(std::ofstream &outfile, const int size0, std::vector<T> &variable) {

    for (int i = 0; i < size0; i++) { outfile << std::setw(15) << variable[i] << "\n"; }
    return;
}

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

void InputParaGriddata(std::ifstream &infile);

inline void OutputMessage(int iview, int istep) {

    std::cout << " ===== Now Computing ===== " << "\n";
    std::cout << "istep:" << std::setw(10) << istep << "\n";
    std::cout << "time:" << std::setw(15) << real_time << std::setw(10) //
              << "iview:" << std::setw(10) << iview << "\n";
    std::cout << "========================== " << "\n";

    return;
}
