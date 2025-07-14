#include "foundation/utils/MathUtils.hpp"
#include <cmath> // For std::sqrt, std::pow, std::atan2, std::abs, std::cyl_bessel_j, std::cyl_neumann

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

} // namespace bem::foundation::utils::MathUtils
