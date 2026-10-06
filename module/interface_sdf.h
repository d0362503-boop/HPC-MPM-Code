#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "module/cal_mat.h"
#include "module/mesh.h"

namespace interface_geometry {

class InterfaceSDF {
  public:
    /**
     * @brief Index the reconstructed surface and finite fluid domains in physical space.
     * @param surface Oriented current solid boundary triangles.
     * @param fluid_domains Current volume-equivalent, grid-aligned fluid particle boxes.
     * @note Contact activates within 1.5 current background-grid spacings.
     */
    InterfaceSDF(std::vector<std::array<std::array<double, 3>, 3>> surface,
                 std::vector<std::array<std::array<double, 3>, 2>> fluid_domains)
        : surface_(std::move(surface)), fluid_domains_(std::move(fluid_domains)) {

        const double search_width = *std::max_element(dxy.begin(), dxy.end());

        for (size_t i = 0; i < this->surface_.size(); i++) {
            std::array<std::array<double, 3>, 2> bounds{this->surface_[i][0], this->surface_[i][0]};
            for (const std::array<double, 3> &vertex : this->surface_[i]) {
                for (int d = 0; d < 3; d++) {
                    bounds[0][d] = std::min(bounds[0][d], vertex[d]);
                    bounds[1][d] = std::max(bounds[1][d], vertex[d]);
                }
            }

            this->IndexBounds(bounds, {search_width, search_width, search_width}, i, this->surface_bins_);
        }

        const std::array<double, 3> contact_spacing{1.5 * dxy[0], 1.5 * dxy[1], 1.5 * dxy[2]};
        for (size_t i = 0; i < this->fluid_domains_.size(); i++) {
            this->IndexBounds(this->fluid_domains_[i], contact_spacing, i, this->fluid_bins_);
        }
    }

    /**
     * @brief Evaluate wet-interface area density chi_w * delta_epsilon(d_s).
     * @param point Current volume quadrature position.
     * @return Interface area per volume [1/m], zero outside the band or on dry surfaces.
     */
    double AreaDensity(const std::array<double, 3> &point) const {

        std::array<double, 3> surface_point;
        std::array<double, 3> normal;
        double signed_distance;
        if (!this->FindSurfacePoint(point, surface_point, signed_distance, normal)) return 0.0;

        std::array<double, 3> normal_spacing;
        for (int d = 0; d < 3; d++) { normal_spacing[d] = normal[d] * dxy[d]; }

        const double band_width = NormVec3(normal_spacing);
        if (std::abs(signed_distance) >= band_width) return 0.0;

        if (!this->IsWetSurfacePoint(surface_point)) return 0.0;

        const double pi = std::acos(-1.0);
        return (1.0 + std::cos(pi * signed_distance / band_width)) / (2.0 * band_width);
    }

  private:
    std::vector<std::array<std::array<double, 3>, 3>> surface_;
    std::vector<std::array<std::array<double, 3>, 2>> fluid_domains_;

    std::map<std::array<int, 3>, std::vector<size_t>> surface_bins_, fluid_bins_;

    /**
     * @brief Activate contact within 1.5 background-grid spacings of the solid surface.
     * @param surface_point Closest physical position on the current solid surface.
     * @return Whether a fluid particle domain lies within the grid-scaled 1.5-distance band.
     */
    bool IsWetSurfacePoint(const std::array<double, 3> &surface_point) const {

        const auto fluid_bin = this->fluid_bins_.find(this->CellIndex(surface_point));
        if (fluid_bin == this->fluid_bins_.end()) return false;

        for (size_t i : fluid_bin->second) {
            std::array<double, 3> separation{};
            for (int d = 0; d < 3; d++) {
                const double nearest =
                    std::clamp(surface_point[d], this->fluid_domains_[i][0][d], this->fluid_domains_[i][1][d]);
                separation[d] = (surface_point[d] - nearest) / (1.5 * dxy[d]);
            }

            // Begin coupling before the finite particle domains touch.
            if (DotVec3(separation, separation) <= 1.0) return true;
        }

        return false;
    }

    /**
     * @brief Locate the closest surface and its signed distance inside the narrow band.
     * @param point Current query position.
     * @param surface_point Closest reconstructed solid surface position.
     * @param signed_distance Distance [m], negative inside the oriented solid surface.
     * @param normal Outward unit normal of the nearest surface triangle.
     * @return Whether the point lies within the reconstructed surface band.
     * @note The sign follows the nearest triangle; area density is even in distance.
     */
    bool FindSurfacePoint(const std::array<double, 3> &point, std::array<double, 3> &surface_point,
                          double &signed_distance, std::array<double, 3> &normal) const {

        const auto bin = this->surface_bins_.find(this->CellIndex(point));
        if (bin == this->surface_bins_.end()) return false;

        const double search_width = *std::max_element(dxy.begin(), dxy.end());
        double closest_distance = search_width;
        for (size_t i : bin->second) {
            const std::array<double, 3> candidate = this->ClosestPointOnTriangle(point, this->surface_[i]);
            const std::array<double, 3> separation = DifferenceVec3(point, candidate);
            const double candidate_distance = NormVec3(separation);

            if (candidate_distance < closest_distance) {
                closest_distance = candidate_distance;
                surface_point = candidate;

                const std::array<std::array<double, 3>, 3> &triangle = this->surface_[i];
                const std::array<double, 3> outward =
                    CrossVec3(DifferenceVec3(triangle[1], triangle[0]), DifferenceVec3(triangle[2], triangle[0]));
                const double normal_length = NormVec3(outward);
                for (int d = 0; d < 3; d++) { normal[d] = outward[d] / normal_length; }

                signed_distance = std::copysign(closest_distance, DotVec3(separation, outward));
            }
        }

        return closest_distance < search_width;
    }

    /**
     * @brief Find the nearest point on a finite triangle, including its edges.
     * @param point Query position in the current configuration.
     * @param triangle Nondegenerate interface triangle.
     * @return Closest physical interface position.
     */
    std::array<double, 3> ClosestPointOnTriangle(const std::array<double, 3> &point,
                                              const std::array<std::array<double, 3>, 3> &triangle) const {

        const std::array<double, 3> ab = DifferenceVec3(triangle[1], triangle[0]);
        const std::array<double, 3> ac = DifferenceVec3(triangle[2], triangle[0]);
        const std::array<double, 3> ap = DifferenceVec3(point, triangle[0]);

        const double ab_squared = DotVec3(ab, ab);
        const double ac_squared = DotVec3(ac, ac);
        const double ab_dot_ac = DotVec3(ab, ac);
        const double ap_dot_ab = DotVec3(ap, ab);
        const double ap_dot_ac = DotVec3(ap, ac);

        // Project onto the triangle plane in barycentric coordinates.
        const std::array<double, 3> normal = CrossVec3(ab, ac);
        const double determinant = DotVec3(normal, normal);
        const double weight_b = (ac_squared * ap_dot_ab - ab_dot_ac * ap_dot_ac) / determinant;
        const double weight_c = (ab_squared * ap_dot_ac - ab_dot_ac * ap_dot_ab) / determinant;

        std::array<double, 3> closest;
        if (weight_b >= 0.0 && weight_c >= 0.0 && weight_b + weight_c <= 1.0) {
            for (int d = 0; d < 3; d++) {
                closest[d] = triangle[0][d] + weight_b * ab[d] + weight_c * ac[d];
            }

            return closest;
        }

        // Outside the triangle, choose the closest point on its three edges.
        double minimum_squared_distance = INFINITY;
        for (int edge = 0; edge < 3; edge++) {
            const std::array<double, 3> &edge_start = triangle[edge];
            const std::array<double, 3> edge_direction = DifferenceVec3(triangle[(edge + 1) % 3], edge_start);
            const double edge_fraction = std::clamp(
                DotVec3(DifferenceVec3(point, edge_start), edge_direction) / DotVec3(edge_direction, edge_direction), 0.0, 1.0);

            std::array<double, 3> candidate;
            for (int d = 0; d < 3; d++) { candidate[d] = edge_start[d] + edge_fraction * edge_direction[d]; }

            const std::array<double, 3> separation = DifferenceVec3(point, candidate);
            const double squared_distance = DotVec3(separation, separation);

            if (squared_distance < minimum_squared_distance) {
                minimum_squared_distance = squared_distance;
                closest = candidate;
            }
        }

        return closest;
    }

    /**
     * @brief Locate an indexing cell without assuming the domain origin.
     * @param point Current physical position.
     * @return Integer spatial-bin coordinates.
     */
    std::array<int, 3> CellIndex(const std::array<double, 3> &point) const {

        return {int(std::floor(point[0] / dxy[0])), int(std::floor(point[1] / dxy[1])),
                int(std::floor(point[2] / dxy[2]))};
    }

    /**
     * @brief Index geometry by its expanded physical bounds.
     * @param bounds Lower and upper physical coordinates.
     * @param margin Search extension in each physical direction.
     * @param index Surface or fluid-domain entry being indexed.
     * @param bins Spatial bins receiving that entry.
     */
    void IndexBounds(std::array<std::array<double, 3>, 2> bounds, const std::array<double, 3> &margin,
                     size_t index, std::map<std::array<int, 3>, std::vector<size_t>> &bins) {

        for (int d = 0; d < 3; d++) {
            bounds[0][d] -= margin[d];
            bounds[1][d] += margin[d];
        }

        const std::array<int, 3> lower = this->CellIndex(bounds[0]);
        const std::array<int, 3> upper = this->CellIndex(bounds[1]);

        for (int k = lower[2]; k <= upper[2]; k++) {
            for (int j = lower[1]; j <= upper[1]; j++) {
                for (int i = lower[0]; i <= upper[0]; i++) { bins[{i, j, k}].push_back(index); }
            }
        }

        return;
    }
};

} // namespace interface_geometry
