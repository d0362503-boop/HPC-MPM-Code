#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace interface_geometry {

class InterfaceSDF {
  public:
    /**
     * @brief Index an oriented triangulated surface for narrow-band distance queries.
     * @param surface Nondegenerate boundary triangles oriented toward the exterior.
     * @param bin_spacing Positive physical bin widths in the three coordinate directions.
     * @note The search band is the largest bin width; no solver or fluid state is required.
     */
    InterfaceSDF(std::vector<std::array<std::array<double, 3>, 3>> surface,
                 const std::array<double, 3> &bin_spacing);

    /**
     * @brief Locate the closest surface and its signed distance inside the narrow band.
     * @param point Current query position.
     * @param surface_point Closest position on the indexed surface.
     * @param signed_distance Distance [m], negative on the nearest triangle's interior side.
     * @param normal Outward unit normal of the nearest surface triangle.
     * @return Whether the point lies within the reconstructed surface band.
     * @note The sign follows the nearest triangle, not a global inside/outside winding test.
     */
    bool FindSurfacePoint(const std::array<double, 3> &point, std::array<double, 3> &surface_point,
                          double &signed_distance, std::array<double, 3> &normal) const;

  private:
    std::vector<std::array<std::array<double, 3>, 3>> surface_;
    std::array<double, 3> bin_spacing_;
    std::map<std::array<int, 3>, std::vector<std::size_t>> surface_bins_;

    /**
     * @brief Find the nearest point on a finite triangle, including its edges.
     * @param point Query position in the current configuration.
     * @param triangle Nondegenerate interface triangle.
     * @return Closest physical interface position.
     */
    std::array<double, 3> ClosestPointOnTriangle(const std::array<double, 3> &point,
                                              const std::array<std::array<double, 3>, 3> &triangle) const;

    /**
     * @brief Locate an indexing cell without assuming the domain origin.
     * @param point Current physical position.
     * @return Integer spatial-bin coordinates.
     */
    std::array<int, 3> CellIndex(const std::array<double, 3> &point) const;

    /**
     * @brief Index geometry by its expanded physical bounds.
     * @param bounds Lower and upper physical coordinates.
     * @param margin Search extension in each physical direction.
     * @param index Surface triangle being indexed.
     */
    void IndexBounds(std::array<std::array<double, 3>, 2> bounds, const std::array<double, 3> &margin,
                     std::size_t index);
};

} // namespace interface_geometry
