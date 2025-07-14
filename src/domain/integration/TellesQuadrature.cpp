#include "domain/integration/TellesQuadrature.hpp"
#include <numbers>

namespace bem::domain::integration {

TellesQuadrature::TellesQuadrature(int order) : order_(order) {
  setupStandardPoints();
}

void TellesQuadrature::setupStandardPoints() {
  standard_points_.clear();
  switch (order_) {
  case 2: {
    Real xi = std::numbers::inv_sqrt3;
    standard_points_.emplace_back(-xi, 1.0, Point2D(-xi, 0.0));
    standard_points_.emplace_back(xi, 1.0, Point2D(xi, 0.0));
    break;
  }
  case 3: {
    Real xi1 = std::sqrt(3.0 / 5.0);
    standard_points_.emplace_back(-xi1, 5.0 / 9.0, Point2D(-xi1, 0.0));
    standard_points_.emplace_back(0.0, 8.0 / 9.0, Point2D(0.0, 0.0));
    standard_points_.emplace_back(xi1, 5.0 / 9.0, Point2D(xi1, 0.0));
    break;
  }
  default:
    throw std::runtime_error("Unsupported Telles quadrature order");
  }
}

Real TellesQuadrature::tellesTransform(Real xi, Real xi_star) {
  return (xi * xi * (3.0 - 2.0 * std::abs(xi_star)) * xi_star / 2.0) + xi;
}

Real TellesQuadrature::tellesTransformDerivative(Real xi, Real xi_star) {
  return 1.0 + (3.0 * xi * (1.0 - std::abs(xi_star)));
}

Real TellesQuadrature::computeXiStar(const Element &element, const Point2D &x) {
  Real dx = element.end.x - element.start.x;
  Real dy = element.end.y - element.start.y;
  Real length2 = (dx * dx) + (dy * dy);
  Real proj =
      ((x.x - element.start.x) * dx + (x.y - element.start.y) * dy) / length2;
  return (2.0 * proj) - 1.0;
}

Complex TellesQuadrature::integrate(
    const Element &element,
    const std::function<Complex(const Point2D &)> &integrand) const {
  Complex result(0.0, 0.0);
  Real xi_star = computeXiStar(element, element.midpoint());
  const double jacobian = element.jacobian();

  for (const auto &qp : standard_points_) {
    Real phi = tellesTransform(qp.xi, xi_star);
    Real dphi_dxi = tellesTransformDerivative(qp.xi, xi_star);
    Point2D y = element.parametricPoint(phi);
    Real weight = qp.weight * dphi_dxi;
    result += weight * integrand(y) * jacobian;
  }
  return result;
}

} // namespace bem::domain::integration
