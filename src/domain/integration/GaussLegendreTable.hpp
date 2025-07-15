#pragma once

#include "types/MathObjectTypes.hpp"

namespace bem::domain::integration {

using types::Complex;
using types::Point2D;
using types::Real;

struct QuadraturePoint {
  Real xi{0.0};     // Parameter coordinate [-1, 1]
  Real weight{0.0}; // Integration weight
  Point2D point;    // Physical coordinate

  QuadraturePoint(Real xi_val, Real weight_val, const Point2D &point_val)
      : xi(xi_val), weight(weight_val), point(point_val) {
  }
};

class GaussLegendreTable {
public:
  static std::vector<QuadraturePoint> getPoints(int order);
};

} // namespace bem::domain::integration
