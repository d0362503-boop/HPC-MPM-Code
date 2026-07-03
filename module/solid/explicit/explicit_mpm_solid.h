#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../solver/crsmat.h"
#include "../solid_material_point.h"

namespace explicitmpm {

class ExplicitSolidMPM : public SolidMaterialPointBase {
  private:
    /**
     * @brief Incremental deformation gradient from the current step, stored so
     *        the position/stress update phase can reuse it.
     */
    std::vector<std::array<std::array<double, 3>, 3>> delta_def_grad_;

  public:
    /**
     * @brief Read and initialize the standalone explicit-solid case input data.
     *
     * Loads global control parameters, mesh/time derived quantities, boundary
     * conditions, and particle data for the explicit solid solver path.
     */
    void DataInput();
    /**
     * @brief Map solid particle mass/momentum/force to control points (P2G).
     */
    void Particle2Node() override;

    /**
     * @brief MUSL velocity projection: re-map particle momentum to control points
     *        and recompute nodal velocity after the G2P update.
     */
    void DoMUSL();

    /**
     * @brief Map updated nodal kinematics back to solid particles (G2P).
     */
    void Node2Particle() override;

  private:
    /**
     * @brief G2P velocity update followed by MUSL projection.
     */
    void G2PVelocityAndMUSL();

    /**
     * @brief Update the particle deformation gradient for the current step.
     *
     * Computes the nodal displacement increment and the incremental deformation
     * gradient, and applies the optional F-bar volumetric correction.
     */
    void UpdateDeformationGradient();

    /**
     * @brief Update particle position, volume, and stress after the deformation
     *        gradient step.
     */
    void UpdateParticlePositionAndStress();

  public:
    /**
     * @brief Compute the F-bar corrected deformation gradient for each particle.
     *
     * Implements the nodal F-bar projection described in:
     *   "Circumventing volumetric locking in explicit material point methods:
     *    A simple, efficient, and general approach"
     *
     * @param delta_def_grad      Incremental deformation-gradient tensor for each particle (modified in-place).
     * @param det_delta_def_grad  Determinant of the incremental deformation gradient for each particle.
     */
    void ComputeDefGradBar(std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                           const std::vector<double> &det_delta_def_grad);

    /**
     * @brief Compute nodal acceleration from nodal force/mass and apply boundary conditions.
     */
    void SolveSolid() override;

};
} // namespace explicitmpm
