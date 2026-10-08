#include "interface_sdf.h"

#include <algorithm>
#include <cmath>
#include <utility>

#include "module/cal_mat.h"

namespace interface_geometry {

InterfaceSDF::InterfaceSDF(std::vector<std::array<std::array<double, 3>, 3>> surface,
                           const std::array<double, 3> &bin_spacing)
    : surface_(std::move(surface)), bin_spacing_(bin_spacing) {

    const double search_width = *std::max_element(this->bin_spacing_.begin(), this->bin_spacing_.end());

    for (std::size_t i = 0; i < this->surface_.size(); i++) {
        std::array<std::array<double, 3>, 2> bounds{this->surface_[i][0], this->surface_[i][0]};
        for (const std::array<double, 3> &vertex : this->surface_[i]) {
            for (int d = 0; d < 3; d++) {
                bounds[0][d] = std::min(bounds[0][d], vertex[d]);
                bounds[1][d] = std::max(bounds[1][d], vertex[d]);
            }
        }

        this->IndexBounds(bounds, {search_width, search_width, search_width}, i);
    }
}

bool InterfaceSDF::FindSurfacePoint(const std::array<double, 3> &point, std::array<double, 3> &surface_point,
                                    double &signed_distance, std::array<double, 3> &normal) const {

    const auto bin = this->surface_bins_.find(this->CellIndex(point));
    if (bin == this->surface_bins_.end()) return false;

    const double search_width = *std::max_element(this->bin_spacing_.begin(), this->bin_spacing_.end());
    double closest_distance = search_width;
    for (std::size_t i : bin->second) {
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

std::array<double, 3> InterfaceSDF::ClosestPointOnTriangle(const std::array<double, 3> &point,
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

std::array<int, 3> InterfaceSDF::CellIndex(const std::array<double, 3> &point) const {

    return {int(std::floor(point[0] / this->bin_spacing_[0])), int(std::floor(point[1] / this->bin_spacing_[1])),
            int(std::floor(point[2] / this->bin_spacing_[2]))};
}

void InterfaceSDF::IndexBounds(std::array<std::array<double, 3>, 2> bounds, const std::array<double, 3> &margin,
                               std::size_t index) {

    for (int d = 0; d < 3; d++) {
        bounds[0][d] -= margin[d];
        bounds[1][d] += margin[d];
    }

    const std::array<int, 3> lower = this->CellIndex(bounds[0]);
    const std::array<int, 3> upper = this->CellIndex(bounds[1]);

    for (int k = lower[2]; k <= upper[2]; k++) {
        for (int j = lower[1]; j <= upper[1]; j++) {
            for (int i = lower[0]; i <= upper[0]; i++) { this->surface_bins_[{i, j, k}].push_back(index); }
        }
    }

    return;
}

} // namespace interface_geometry
