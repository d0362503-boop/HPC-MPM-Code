#pragma once

#include <fstream>
#include <string>
#include <vector>

struct Mesh;
struct ControlMesh;
struct Grid;

class BoundaryCondition {
  public:
    int ibc = 0;
    std::vector<int> nbc;
    std::vector<double> fbc;
    std::vector<int> global_nbc;    // immutable global IDs
    std::vector<double> global_fbc; // immutable prescribed values

    /**
     * @brief Read boundary-condition node IDs and prescribed values from input stream.
     * @param infile  Input file stream.
     * @param if_fbc  If true, also read prescribed values into `fbc`.
     */
    void BCInput(std::ifstream &infile, bool if_fbc = true);

    /**
     * @brief Write boundary-condition node IDs and values to an output stream.
     * @param outfile   Output file stream.
     * @param bc_name   Name of the boundary condition (used in the header).
     * @param if_fbc    If true, also write prescribed values.
     */
    void BCOutput(std::ofstream &outfile, std::string bc_name, bool if_fbc = true);

    /**
     * @brief Apply prescribed boundary values `scale * fbc` to `variable` at constrained DOFs.
     * @param nn       Offset into `variable` for the component being constrained.
     * @param variable Vector to be modified in-place.
     * @param scale    Load factor applied to the prescribed values (the current `facl`).
     */
    void BCSetVal(int nn, std::vector<double> &variable, double scale);

    /**
     * @brief Zero out `variable` at constrained DOFs.
     * @param nn       Offset into `variable` for the component being constrained.
     * @param variable Vector to be modified in-place.
     */
    void BCSetZero(int nn, std::vector<double> &variable);

    /**
     * @brief Set the displacement increment `scale * fbc` at constrained DOFs.
     * @param nn       Offset into `variable` for the component being constrained.
     * @param variable Displacement increment modified in-place.
     * @param scale    Velocity-to-increment factor `facl * dt` of the current step.
     * @note Matches the Newmark endpoint velocity `facl * fbc` because the step-start
     * `nvel` and `naccel` already hold `facl * fbc` and zero there (ApplyVelocityBC / ApplyAccelerationBC in P2G).
     */
    void BCSetDispl(int nn, std::vector<double> &variable, double scale);

    /**
     * @brief Remove the wall-normal component of a nodal vector at constrained control points (free slip).
     * @param variable Nodal vector field (x, y, z blocks) modified in-place.
     * @param normal   Nodal unit normal field in the same block layout.
     * @param cmesh    Control points of this rank (component offsets are used).
     * @param tol      Normal magnitude below which the whole vector is zeroed instead.
     */
    void BCSetFreeSlip(std::vector<double> &variable, const std::vector<double> &normal, const ControlMesh &cmesh,
                       double tol);

    /**
     * @brief Cache this rank's control-point BCs as a global MPI-wide BC list.
     *
     * Local control-point IDs are converted to global IDs and gathered from all
     * ranks.  The resulting data remains valid after a DLB repartition.
     * @param grid Grid of this rank (element window and control-point sizes).
     */
    void CaptureGlobalControlPointBC(const Grid &grid);

    /**
     * @brief Rebuild local control-point BC IDs for the current MPI region.
     * @param grid Grid of this rank after the repartition.
     */
    void RebuildLocalControlPointBC(const Grid &grid);

    /**
     * @brief Cache this rank's element BCs as a global MPI-wide BC list.
     *
     * This is used by fluid inflow BCs, whose IDs refer to background elements.
     * @param mesh Background mesh of this rank (element window and counts).
     */
    void CaptureGlobalElementBC(const Mesh &mesh);

    /**
     * @brief Rebuild local background-element BC IDs for the current MPI region.
     * @param mesh Background mesh of this rank after the repartition.
     */
    void RebuildLocalElementBC(const Mesh &mesh);

  private:
    /**
     * @brief Gather and de-duplicate global BC entries from all MPI ranks.
     * @param local_global_ids Global IDs represented by this rank's local BC entries.
     * @param has_values Whether each BC entry has a prescribed scalar value.
     */
    void CacheGlobalEntries(const std::vector<int> &local_global_ids, bool has_values);
};
