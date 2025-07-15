#include "GaussQuadrature.hpp"
#include "GaussLegendreTable.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include <cmath>
#include <vector>

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Element;
using bem::types::Point2D;

GaussQuadrature::GaussQuadrature(int order) : order_(order) {
  setupQuadraturePoints();
}

Complex GaussQuadrature::integrate(
    const Element &element,
    const Point2D & /* x */,
    const std::function<Complex(const Point2D &)> &integrand) const {
  Complex result(0.0, 0.0);
  const double length = element.length();
  for (const auto &point : points_) {
    Point2D y = element.parametricPoint(point.xi);
    result += point.weight * integrand(y) * (length / 2);
  }
  return result;
}

void GaussQuadrature::setupQuadraturePoints() {
  points_ = GaussLegendreTable::getPoints(order_);
}

} // namespace bem::domain::integration
