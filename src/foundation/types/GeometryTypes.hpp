#pragma once

#include "ComplexTypes.hpp"

namespace bem::types {

// Boundary element types
enum class ElementType : std::uint8_t {
  LINEAR,    // 2-node linear element
  QUADRATIC, // 3-node quadratic element
  CUBIC      // 4-node cubic element
};

// Boundary condition types
enum class BoundaryConditionType : std::uint8_t {
  DIRICHLET, // Specified potential
  NEUMANN,   // Specified normal derivative
  ROBIN,     // Mixed boundary condition
  IMPEDANCE  // Acoustic impedance boundary
};

// Integration rule types
enum class IntegrationType : std::uint8_t {
  GAUSS,     // Gauss-Legendre quadrature
  ADAPTIVE,  // Adaptive quadrature
  ANALYTICAL // Analytical integration (when possible)
};

// Line segment for boundary elements
struct LineSegment {
  Point2D start;
  Point2D end;

  LineSegment() = default;
  LineSegment(const Point2D &p1, const Point2D &p2) : start(p1), end(p2) {}

  [[nodiscard]] Real length() const;
  [[nodiscard]] Point2D midpoint() const;
  [[nodiscard]] Vector2D tangent() const;
  [[nodiscard]] Vector2D normal() const;                // Outward normal
  [[nodiscard]] Point2D parametricPoint(Real xi) const; // xi ∈ [-1, 1]
  [[nodiscard]] Real jacobian() const;                  // For integration
};

// Quadrature point
struct QuadraturePoint {
  Real xi{0.0};     // Parameter coordinate [-1, 1]
  Real weight{0.0}; // Integration weight
  Point2D point;    // Physical coordinate

  QuadraturePoint() = default;
  QuadraturePoint(Real xi_val, Real weight_val, const Point2D &point_val)
      : xi(xi_val), weight(weight_val), point(point_val) {}
};

} // namespace bem::types
