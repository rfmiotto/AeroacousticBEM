#pragma once

#include "MathObjectTypes.hpp"

namespace bem::types {

// Boundary condition types
enum class BoundaryConditionType : std::uint8_t {
  DIRICHLET, // Specified potential
  NEUMANN,   // Specified normal derivative
  ROBIN,     // Mixed boundary condition
  IMPEDANCE, // Acoustic impedance boundary
};

// Boundary value
struct BoundaryValue {
  Complex potential;         // φ (acoustic potential)
  Complex normal_derivative; // ∂φ/∂n (normal velocity)
  bool potential_known{false};
  bool derivative_known{false};

  BoundaryValue() = default;
  explicit BoundaryValue(Complex pot, bool pot_known = true)
      : potential(pot), potential_known(pot_known) {
  }
};

// Function types for boundary conditions
using BoundaryFunction = std::function<Complex(const Point2D &)>;
using NormalFunction = std::function<Vector2D(const Point2D &)>;

} // namespace bem::types
