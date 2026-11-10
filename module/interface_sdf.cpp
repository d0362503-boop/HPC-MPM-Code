#include "interface_sdf.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "module/cal_mat.h"
#include "module/mesh.h"
#include "module/shape_function.h"

namespace interface_geometry {

InterfaceSDF::InterfaceSDF(std::vector<double> solid_phi) : solid_phi_(std::move(solid_phi)) {}

double InterfaceSDF::Evaluate(const std::array<double, 3> &point, std::array<double, 3> &gradient) const {
    gradient = {};
    for (double coordinate : point) {
        if (!std::isfinite(coordinate)) { return std::numeric_limits<double>::quiet_NaN(); }
    }

    int cell;
    if (!LocateLocalElement(point, cell)) return 0.0;

    int count;
    std::vector<int> nodes;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;
    MakeSF(cell, point, idimc, xynodec, nodes, count, sf, dsf);

    // Subtract the target before interpolation to preserve constant levels.
    double level = 0.0;
    for (int i = 0; i < count; ++i) {
        const double value = this->solid_phi_[nodes[i]] - 0.5;
        level += sf[i] * value;
        for (int d = 0; d < 3; ++d) { gradient[d] += dsf[i][d] * value; }
    }
    return 0.5 + level;
}

bool InterfaceSDF::FindSurfacePoint(const std::array<double, 3> &point,
                                    std::array<double, 3> &surface_point) const {
    std::array<double, 3> gradient;
    const double level = this->Evaluate(point, gradient) - 0.5;
    if (!(level > 0.0)) return false;

    const double gradient_squared = DotVec3(gradient, gradient);
    const double min_spacing = *std::min_element(dxy.begin(), dxy.end());
    const double min_gradient = 1.0e-12 / min_spacing;
    if (!std::isfinite(gradient_squared) || gradient_squared <= min_gradient * min_gradient) return false;

    const double scale = level / gradient_squared;
    for (int d = 0; d < 3; ++d) {
        surface_point[d] = std::clamp(point[d] - scale * gradient[d], xyminw[d], xymaxw[d]);
    }
    return true;
}

} // namespace interface_geometry
