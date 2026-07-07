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
  public:
    /**
     * @brief Read and initialize standalone explicit-solid input data.
     */
    void DataInput();

    /**
     * @brief Map solid particle mass/volume/momentum/force to control points (P2G).
     */
    void Particle2Node() override;

    /**
     * @brief Map updated nodal kinematics back to solid particles (G2P).
     */
    void Node2Particle() override;

    /**
     * @brief Compute nodal acceleration from force/mass and apply BCs.
     */
    void SolveSolid() override;

  private:
    // Incremental deformation gradient for the current step.
    std::vector<std::array<std::array<double, 3>, 3>> delta_def_grad;

    /**
     * @brief G2P velocity update followed by MUSL projection.
     */
    void G2PVelocity();

    /**
     * @brief MUSL velocity projection.
     */
    void DoMUSL();

    /**
     * @brief Update particle deformation gradient.
     */
    void UpdateDeformationGradient();

    /**
     * @brief Update particle position, volume, and stress.
     */
    void UpdateParticlePositionAndStress();

    /**
     * @brief F-bar corrected deformation gradient.
     *
     * Nodal F-bar projection from:
     *   "Circumventing volumetric locking in explicit material point methods".
     *
     * @param delta_def_grad      Incremental deformation gradient (modified in-place).
     * @param det_delta_def_grad  Determinant of incremental deformation gradient.
     */
    void ComputeDefGradBar(const std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                           const std::vector<double> &det_delta_def_grad);
};

} // namespace explicitmpm
