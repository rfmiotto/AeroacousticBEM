#include "foundation/utils/MathUtils.hpp"
#include <cmath> // For std::sqrt, std::pow, std::atan2, std::abs, std::cyl_bessel_j, std::cyl_neumann
#include <numeric> // For std::inner_product (potentially, or manual dot product)
#include <stdexcept> // For std::invalid_argument
#include <vector>    // For std::vector

namespace bem::foundation::utils::MathUtils {

double distance(const types::Point2D &p1, const types::Point2D &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return std::sqrt((dx * dx) + (dy * dy));
}

double distanceSquared(const types::Point2D &p1, const types::Point2D &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);
}

double angle(const types::Point2D &p1, const types::Point2D &p2) {
  return std::atan2(p2.y - p1.y, p2.x - p1.x);
}

bool isZero(double value, double tol) {
  return std::abs(value) < tol;
}

bool isEqual(double v1, double v2, double tol) {
  return std::abs(v1 - v2) < tol;
}

double magnitude(types::Complex z) {
  return std::abs(z);
}

double phase(types::Complex z) {
  return std::arg(z);
}

// std::pair<std::vector<double>, std::vector<double>>
// gaussLegendreQuadrature(int n_points) {
//   if (n_points == 2) {
//     double point = std::numbers::inv_sqrt3;
//     std::vector<double> points = {-point, point};
//     std::vector<double> weights = {1.0, 1.0};
//     return {points, weights};
//   } else {
//     throw std::invalid_argument("gaussLegendreQuadrature is currently only "
//                                 "implemented for n_points = 2.");
//   }
// }
//
// std::pair<types::Point2D, double>
// lineIntegrationTransform(const types::Point2D &start, const types::Point2D
// &end,
//                          double xi) {
//   double common_factor_start = (1.0 - xi) / 2.0;
//   double common_factor_end = (1.0 + xi) / 2.0;
//
//   types::Point2D transformed_point;
//   transformed_point.x =
//       common_factor_start * start.x + common_factor_end * end.x;
//   transformed_point.y =
//       common_factor_start * start.y + common_factor_end * end.y;
//
//   double jacobian = distance(start, end) / 2.0;
//
//   return {transformed_point, jacobian};
// }
//
} // namespace bem::foundation::utils::MathUtils
