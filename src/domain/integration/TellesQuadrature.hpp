#pragma once

#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include <cmath>
#include <vector>

namespace bem::domain::integration {

using bem::types::Panel;
using bem::types::Point2D;
using bem::types::QuadraturePoint;
using bem::types::Real;

/**
 * @brief Telles quadrature for weakly singular integrals
 *        (see: Telles, 1987 - A self-adaptive co-ordinate transformation)
 */
class TellesQuadrature {
public:
  explicit TellesQuadrature(int order) : order_(order) {
    setupStandardPoints();
  }

  /**
   * @brief Compute transformed quadrature points using Telles' method
   * @param panel Panel over which integration is done
   * @param x Field point where the singularity is located (typically on the
   * panel)
   * @return Vector of transformed QuadraturePoint objects
   */
  [[nodiscard]] std::vector<QuadraturePoint> compute(const Panel &panel,
                                                     const Point2D &x) const {
    std::vector<QuadraturePoint> transformed;
    transformed.reserve(points_.size());

    // Compute projection of x onto panel's local parametric coordinate (xi*)
    const Real xi_star = computeXiStar(panel, x);

    // Telles transformation
    for (const auto &qp : points_) {
      const Real xi = qp.xi;
      Real phi = tellesTransform(xi, xi_star);
      Real dphi_dxi = tellesTransformDerivative(xi, xi_star);

      Point2D y = panel.parametricPoint(phi);
      Real weight = qp.weight * dphi_dxi;

      transformed.emplace_back(phi, weight, y);
    }

    return transformed;
  }

private:
  int order_;
  std::vector<QuadraturePoint> points_;

  void setupStandardPoints(); // Gauss-Legendre in [-1,1]

  // Transformations
  [[nodiscard]] static Real tellesTransform(Real xi, Real xi_star) {
    return (xi * xi * (3 - 2 * std::abs(xi_star)) * xi_star / 2.0) + xi;
  }

  [[nodiscard]] static Real tellesTransformDerivative(Real xi, Real xi_star) {
    return 1.0 + (3.0 * xi * (1.0 - std::abs(xi_star)));
  }

  [[nodiscard]] static Real computeXiStar(const Panel &panel,
                                          const Point2D &x) {
    // Solve (1 - xi)/2 * start + (1 + xi)/2 * end = x => solve for xi
    Real dx = panel.end.x - panel.start.x;
    Real dy = panel.end.y - panel.start.y;
    Real length2 = (dx * dx) + (dy * dy);
    Real proj =
        ((x.x - panel.start.x) * dx + (x.y - panel.start.y) * dy) / length2;
    return (2.0 * proj) - 1.0; // map to [-1,1]
  }
};

inline void TellesQuadrature::setupStandardPoints() {
  points_.clear();
  switch (order_) {
  case 2: {
    Real xi = std::numbers::inv_sqrt3;
    points_.emplace_back(-xi, 1.0, Point2D(-xi, 0.0));
    points_.emplace_back(xi, 1.0, Point2D(xi, 0.0));
    break;
  }
  case 3: {
    Real xi1 = std::sqrt(3.0 / 5.0);
    points_.emplace_back(-xi1, 5.0 / 9.0, Point2D(-xi1, 0.0));
    points_.emplace_back(0.0, 8.0 / 9.0, Point2D(0.0, 0.0));
    points_.emplace_back(xi1, 5.0 / 9.0, Point2D(xi1, 0.0));
    break;
  }
  default:
    throw std::runtime_error("Unsupported Telles quadrature order");
  }
}

} // namespace bem::domain::integration
