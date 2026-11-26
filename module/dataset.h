#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

inline int ista, iend, iout, istep, nlstep, rstflag, nmat, npropmax;
inline double dt, dti, real_time, dlstep, facl = 1.0e0;
inline double mtol, dtol;
inline std::string parafile, gridfile, outfile, pointfile;
inline std::array<int, 3> npxye{};
inline std::array<double, 3> bb{};
inline std::vector<double> matflag, iprop, nprop;
inline std::vector<std::vector<double>> mat_prop;
inline constexpr std::array<std::array<int, 3>, 3> idm = {{{{0, 5, 4}}, {{5, 1, 3}}, {{4, 3, 2}}}}; // Voigt index map
inline constexpr std::array<std::array<double, 3>, 3> eye_mat = {                                // 3x3 identity
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
 * @param init_value Value assigned to each entry after resizing; defaults to zero for numbers and std::array.
 */
template <typename T> void VectorAssign(int size, std::vector<T> &variable, T init_value = T()) {

    variable.assign(size, init_value);

    return;
}

/**
 * @brief Initialize PETSc and the global MPI rank information.
 * @param argc Number of command-line arguments, updated by PETSc.
 * @param argv Command-line arguments containing runtime solver options.
 * @return Wall-clock start time in seconds for this simulation.
 */
double InitializeSimulation(int &argc, char **&argv);

/**
 * @brief Report elapsed simulation time and finalize PETSc.
 * @param start_time Wall-clock start time returned by InitializeSimulation.
 */
void FinalizeSimulation(double start_time);
