#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

inline int ista, iend, iout, istep, nlstep, rstflag, nmat, npropmax;
inline double dt, dti, real_time, dlstep, facl = 1.0e0;
inline double mtol, dtol;
inline std::string parafile, gridfile, outfile, pointfile;
inline std::vector<int> npxye(3);
inline std::vector<double> bb(3), matflag, iprop, nprop;
inline std::vector<std::vector<double>> mat_prop;
inline std::array<std::array<int, 3>, 3> idm = {{{{0, 5, 4}}, {{5, 1, 3}}, {{4, 3, 2}}}};
inline std::array<std::array<double, 3>, 3> eye_mat = {
    {{{1.0e0, 0.0e0, 0.0e0}}, {{0.0e0, 1.0e0, 0.0e0}}, {{0.0e0, 0.0e0, 1.0e0}}}};

enum class MapScheme { PIC, FLIP, TPIC, APIC };

inline MapScheme ParseMapScheme(const std::string &s) {
    if (s == "FLIP") return MapScheme::FLIP;
    if (s == "PIC") return MapScheme::PIC;
    if (s == "TPIC") return MapScheme::TPIC;
    if (s == "APIC") return MapScheme::APIC;
    throw std::invalid_argument("ParseMapScheme: unknown mapping scheme '" + s + "'");
}

/**
 * @brief Resize a vector and fill every entry with the given initial value.
 * @param size Target number of entries.
 * @param variable Vector to reset.
 * @param init_value Value assigned to each entry after resizing.
 */
template <typename T> void VectorAssign(int size, std::vector<T> &variable, T init_value = T()) {

    variable.assign(size, init_value);

    return;
}

/**
 * @brief Resize a vector of fixed-size arrays and value-initialize every entry.
 * @param size Target number of entries.
 * @param variable Vector of arrays to reset.
 */
template <typename T, std::size_t N> void VectorAssign(int size, std::vector<std::array<T, N>> &variable) {

    variable.assign(size, std::array<T, N>{});

    return;
}
