#include "integration/GaussQuadrature.hpp"
#include "foundation/exceptions/BEMException.hpp"
#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include <array>
#include <cmath>
#include <vector>

namespace bem::domain::integration {

using bem::foundation::exceptions::BEMException;
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
  points_.clear();
  switch (order_) {
  case 1:
    points_.emplace_back(0.0, 2.0, Point2D(0.0, 0.0));
    break;
  case 2: {
    double xi = std::numbers::inv_sqrt3;
    points_.emplace_back(-xi, 1.0, Point2D(-xi, 0.0));
    points_.emplace_back(xi, 1.0, Point2D(xi, 0.0));
    break;
  }
  case 3: {
    double xi = std::sqrt(3.0 / 5.0);
    points_.emplace_back(-xi, 5.0 / 9.0, Point2D(-xi, 0.0));
    points_.emplace_back(0.0, 8.0 / 9.0, Point2D(0.0, 0.0));
    points_.emplace_back(xi, 5.0 / 9.0, Point2D(xi, 0.0));
    break;
  }
  case 4: {
    double xi1 = std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
    double xi2 = std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
    double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
    double w2 = (18.0 - std::sqrt(30.0)) / 36.0;
    points_.emplace_back(-xi2, w2, Point2D(-xi2, 0.0));
    points_.emplace_back(-xi1, w1, Point2D(-xi1, 0.0));
    points_.emplace_back(xi1, w1, Point2D(xi1, 0.0));
    points_.emplace_back(xi2, w2, Point2D(xi2, 0.0));
    break;
  }
  case 5: {
    double xi1 = 0.0;
    double xi2 = std::sqrt(5.0 - (2.0 * std::sqrt(10.0 / 7.0))) / 3.0;
    double xi3 = std::sqrt(5.0 + (2.0 * std::sqrt(10.0 / 7.0))) / 3.0;
    double w1 = 128.0 / 225.0;
    double w2 = (322.0 + 13.0 * std::sqrt(70.0)) / 900.0;
    double w3 = (322.0 - 13.0 * std::sqrt(70.0)) / 900.0;
    points_.emplace_back(-xi3, w3, Point2D(-xi3, 0.0));
    points_.emplace_back(-xi2, w2, Point2D(-xi2, 0.0));
    points_.emplace_back(xi1, w1, Point2D(xi1, 0.0));
    points_.emplace_back(xi2, w2, Point2D(xi2, 0.0));
    points_.emplace_back(xi3, w3, Point2D(xi3, 0.0));
    break;
  }
  case 6: {
    double xi1 = 0.2386191860831969;
    double xi2 = 0.6612093864662645;
    double xi3 = 0.9324695142031521;
    double w1 = 0.4679139345726910;
    double w2 = 0.3607615730481386;
    double w3 = 0.1713244923791703;
    points_.emplace_back(-xi3, w3, Point2D(-xi3, 0.0));
    points_.emplace_back(-xi2, w2, Point2D(-xi2, 0.0));
    points_.emplace_back(-xi1, w1, Point2D(-xi1, 0.0));
    points_.emplace_back(xi1, w1, Point2D(xi1, 0.0));
    points_.emplace_back(xi2, w2, Point2D(xi2, 0.0));
    points_.emplace_back(xi3, w3, Point2D(xi3, 0.0));
    break;
  }
  case 8: {
    std::array<double, 4> xi = {0.1834346424956498, 0.5255324099163290,
                                0.7966664774136267, 0.9602898564975363};
    std::array<double, 4> w = {0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
    for (int i = 0; i < 4; ++i) {
      points_.emplace_back(-xi[i], w[i], Point2D(-xi[i], 0.0));
      points_.emplace_back(xi[i], w[i], Point2D(xi[i], 0.0));
    }
    break;
  }
  default:
    throw BEMException("Unsupported quadrature order: " +
                       std::to_string(order_));
  }
}

} // namespace bem::domain::integration
