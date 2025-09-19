#include "foundation/utils/MathUtils.hpp"
#include "types/MathObjectTypes.hpp"
#include <cmath> // For std::sqrt, std::pow, std::atan2, std::abs, std::cyl_bessel_j, std::cyl_neumann

namespace bem::foundation::utils::MathUtils {

using types::Complex;
using types::Point2D;

double distance(const Point2D &p1, const Point2D &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return std::sqrt((dx * dx) + (dy * dy));
}

double distanceSquared(const Point2D &p1, const Point2D &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);
}

double angle(const Point2D &p1, const Point2D &p2) {
  return std::atan2(p2.y - p1.y, p2.x - p1.x);
}

bool isZero(double value, double tol) {
  return std::abs(value) < tol;
}

bool isEqual(double v1, double v2, double tol) {
  return std::abs(v1 - v2) < tol;
}

bool pointsEqual(const Point2D &p1, const Point2D &p2, double tol) {
  return (std::abs(p1.x - p2.x) <= tol) && (std::abs(p1.y - p2.y) <= tol);
}

double magnitude(Complex z) {
  return std::abs(z);
}

double phase(Complex z) {
  return std::arg(z);
}

} // namespace bem::foundation::utils::MathUtils
