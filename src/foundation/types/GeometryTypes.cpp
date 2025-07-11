#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/MathUtils.hpp"

namespace bem::types {
using namespace bem::foundation::utils::MathUtils;

Real Panel::length() const {
  return distance(start, end);
}

Point2D Panel::midpoint() const {
  return Point2D{0.5 * (start.x + end.x), 0.5 * (start.y + end.y)};
}

Vector2D Panel::tangent() const {
  return Vector2D{end.x - start.x, end.y - start.y};
}

Vector2D Panel::normal() const {
  Vector2D t = tangent();
  return Vector2D{-t.y, t.x};
}

// Parametric point with xi ∈ [-1, 1]
Point2D Panel::parametricPoint(Real xi) const {
  const Real a = 0.5 * (1.0 - xi);
  const Real b = 0.5 * (1.0 + xi);
  return Point2D{(a * start.x) + (b * end.x), (a * start.y) + (b * end.y)};
}

// Jacobian for numerical integration in the [-1, 1] parametric domain
Real Panel::jacobian() const {
  return 0.5 * length();
}

} // namespace bem::types
