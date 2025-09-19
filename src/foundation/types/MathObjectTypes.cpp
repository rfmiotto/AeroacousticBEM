#include "MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"
#include <cmath> // for std::abs

namespace bem::types {

// Arithmetic operators
Point2D Point2D::operator+(const Point2D &other) const {
  return {x + other.x, y + other.y};
}

Point2D Point2D::operator-(const Point2D &other) const {
  return {x - other.x, y - other.y};
}

Point2D Point2D::operator*(Real scalar) const {
  return {x * scalar, y * scalar};
}

Point2D Point2D::operator/(Real scalar) const {
  return {x / scalar, y / scalar};
}

// Comparison operators (tolerant)
bool Point2D::operator==(const Point2D &other) const {
  constexpr Real TOL = bem::foundation::utils::Constants::GEOMETRY_TOLERANCE;
  return std::abs(x - other.x) <= TOL && std::abs(y - other.y) <= TOL;
}

bool Point2D::operator!=(const Point2D &other) const {
  return !(*this == other);
}

} // namespace bem::types
