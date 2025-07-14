#pragma once

#include "IQuadratureRule.hpp"
#include <cmath>
#include <vector>

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Point2D;
using bem::types::Real;

class TellesQuadrature : public IQuadratureRule {
public:
  explicit TellesQuadrature(int order);

  [[nodiscard]]
  Complex integrate(
      const Element &element,
      const Point2D &x,
      const std::function<Complex(const Point2D &)> &integrand) const override;

  [[nodiscard]] const std::vector<QuadraturePoint> &getPoints() const {
    return standard_points_;
  }

private:
  int order_;
  std::vector<QuadraturePoint> standard_points_;

  void setupStandardPoints();
  static Real tellesTransform(Real xi, Real xi_star);
  static Real tellesTransformDerivative(Real xi, Real xi_star);
  static Real computeXiStar(const Element &element, const Point2D &x);
};

} // namespace bem::domain::integration
