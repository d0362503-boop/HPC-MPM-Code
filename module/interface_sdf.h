#pragma once

#include <array>
#include <vector>

namespace interface_geometry {

/**
 * @brief Query and project onto the current solid phi=0.5 interface.
 * @note Uses the rank-local background mesh and synchronized control-point coefficients.
 */
class InterfaceSDF {
  public:
    /**
     * @brief Store the solid phase coefficients in rank-local control-point order.
     * @param solid_phi Final solid volume fractions, indexed with x varying fastest.
     */
    explicit InterfaceSDF(std::vector<double> solid_phi);

    /**
     * @brief Evaluate the original phase field and its physical spatial gradient.
     * @param point Physical query position.
     * @param gradient Receives the phase gradient [1/m].
     * @return Solid volume fraction; zero outside this rank and NaN for nonfinite coordinates.
     */
    double Evaluate(const std::array<double, 3> &point, std::array<double, 3> &gradient) const;

    /**
     * @brief Apply one normal correction toward the original phi=0.5 level set.
     * @param point Final fluid centre before geometric correction.
     * @param surface_point Receives the corrected position clipped to the global mesh.
     * @return Whether an interior centre with a usable phase gradient was corrected.
     * @note A curved interface generally leaves a phase residual after this single update.
     */
    bool FindSurfacePoint(const std::array<double, 3> &point, std::array<double, 3> &surface_point) const;

  private:
    std::vector<double> solid_phi_;
};

} // namespace interface_geometry
