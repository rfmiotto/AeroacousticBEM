#pragma once

#include "MathObjectTypes.hpp"

namespace bem::types {

/**
 * @brief Represents the physical region associated with a boundary segment.
 */
enum class SegmentRegion : std::uint8_t {
  UPSTREAM,
  AIRFOIL,
  DOWNSTREAM,
};

// Line segment for boundary elements
struct LineSegment {
  Point2D start;
  Point2D end;
  SegmentRegion region = SegmentRegion::UPSTREAM; // default initialization

  LineSegment() = default;
  LineSegment(const Point2D &p1, const Point2D &p2) : start(p1), end(p2) {
  }

  [[nodiscard]] Real length() const;
  [[nodiscard]] Point2D midpoint() const;
  [[nodiscard]] Vector2D tangent() const;
  [[nodiscard]] Vector2D normal() const;                // Outward normal
  [[nodiscard]] Point2D parametricPoint(Real xi) const; // xi ∈ [-1, 1]
  [[nodiscard]] Real jacobian() const;                  // For integration
};

} // namespace bem::types
