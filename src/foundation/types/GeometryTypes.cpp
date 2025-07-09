#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/MathUtils.hpp"

namespace bem::types {
using namespace bem::foundation::utils::MathUtils;

Real LineSegment::length() const {
  return distance(start, end);
}

Point2D LineSegment::midpoint() const {
  return Point2D{0.5 * (start.x + end.x), 0.5 * (start.y + end.y)};
}

Vector2D LineSegment::tangent() const {
  return Vector2D{end.x - start.x, end.y - start.y};
}

Vector2D LineSegment::normal() const {
  Vector2D t = tangent();
  return Vector2D{-t.y, t.x};
}

// Parametric point with xi ∈ [-1, 1]
Point2D LineSegment::parametricPoint(Real xi) const {
  const Real a = 0.5 * (1.0 - xi);
  const Real b = 0.5 * (1.0 + xi);
  return Point2D{(a * start.x) + (b * end.x), (a * start.y) + (b * end.y)};
}

// Jacobian for numerical integration in the [-1, 1] parametric domain
Real LineSegment::jacobian() const {
  return 0.5 * length();
}

} // namespace bem::types
