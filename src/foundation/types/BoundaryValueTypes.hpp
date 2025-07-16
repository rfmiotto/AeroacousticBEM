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
  Complex potential{0.0};         // φ
  Complex normal_derivative{0.0}; // ∂φ/∂n
};

// Function types for boundary conditions
using BoundaryFunction = std::function<Complex(const Point2D &)>;
using NormalFunction = std::function<Vector2D(const Point2D &)>;

} // namespace bem::types
