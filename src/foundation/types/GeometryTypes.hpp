#pragma once

#include "MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"

namespace bem::types {

/**
 * @brief Represents the physical region associated with a boundary element.
 */
enum class ElementRegion : std::uint8_t {
  UPSTREAM,
  AIRFOIL,
  DOWNSTREAM,
};

// Boundary element types
enum class ElementType : std::uint8_t {
  CONSTANT,  // 1-node constant element
  LINEAR,    // 2-node linear element
  QUADRATIC, // 3-node quadratic element
  CUBIC      // 4-node cubic element
};

// Line segment for boundary elements
struct Element {
  Point2D start;
  Point2D end;
  ElementRegion region = ElementRegion::UPSTREAM;
  ElementType type = ElementType::CONSTANT;

  Element() = default;
  Element(const Point2D &p1, const Point2D &p2) : start(p1), end(p2) {
  }

  [[nodiscard]] Real length() const;
  [[nodiscard]] Point2D midpoint() const;
  [[nodiscard]] Vector2D tangent() const;
  [[nodiscard]] Vector2D normal() const;                // Outward normal
  [[nodiscard]] Point2D parametricPoint(Real xi) const; // xi ∈ [-1, 1]
  [[nodiscard]] Real jacobian() const;                  // For integration
  [[nodiscard]] bool contains(
      const Point2D &pt,
      Real tol = bem::foundation::utils::Constants::GEOMETRY_TOLERANCE) const;
};

} // namespace bem::types
