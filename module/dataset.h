#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

/// Run parameters read from file.dat and the parameter file; fixed during a run.
struct RunParams {
    int ista = 0, iend = 0, iout = 0;                    // step range and output interval
    int nlstep = 0, rstflag = 0;                         // load-ramp steps, restart mode
    double dt = 0.0, dti = 0.0, dlstep = 0.0;            // time step, 1/dt, 1/nlstep
    double mtol = 0.0, dtol = 0.0;                       // small-mass / distance tolerances
    std::string parafile, gridfile, outfile, pointfile; // input / output prefixes
    std::array<int, 3> npxye{};                         // points per element direction
    std::array<double, 3> bb{};                         // body force
};

/// Time-step state advanced by the driver loop.
struct TimeState {
    int istep = 0;          // current time step
    double real_time = 0.0; // current physical time
    double facl = 1.0e0;    // load ramp factor
};

/// Solid material table read from the parameter file.
struct MaterialTable {
    int nmat = 0, npropmax = 0;                  // material count, max properties
    std::vector<double> matflag, iprop, nprop;   // per-material model ids / counts
    std::vector<std::vector<double>> mat_prop;   // per-material properties
};

inline RunParams g_params; // transitional global instance (still read by CrsMat)
inline TimeState g_time;   // transitional global instance (still read by CrsMat)

// ----- Transitional aliases still read by CrsMat; remove together with g_params/g_time -----
inline double &mtol = g_params.mtol;
inline int &istep = g_time.istep;
// -----------------------------------------------------------------------------------------
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
