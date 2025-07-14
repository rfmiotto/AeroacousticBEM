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

bool Panel::contains(const Point2D &pt, Real tol) const {
  // Panel vector
  const Real dx = end.x - start.x;
  const Real dy = end.y - start.y;

  // Vector from panel start to point
  const Real dxp = pt.x - start.x;
  const Real dyp = pt.y - start.y;

  // If the cross product between these two vectors is not zero,
  // it is not aligned → point is not on the line segment
  const Real cross = (dx * dyp) - (dy * dxp);

  if (std::abs(cross) > tol) {
    return false;
  }

  // Dot produt to check if point is within panel extremes
  const Real dot = (dx * dxp) + (dy * dyp);
  const Real len2 = (dx * dx) + (dy * dy);

  return (dot >= -tol) && (dot <= len2 + tol);
}
} // namespace bem::types
