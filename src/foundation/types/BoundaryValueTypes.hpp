#pragma once

#include "MathObjectTypes.hpp"

namespace bem::types {

// Boundary condition types
enum class BCType : std::uint8_t {
  DIRICHLET, // Specified potential
  NEUMANN,   // Specified normal derivative
  KUTTA,     // Kutta condition
};

// Boundary value
struct BoundaryCondition {
  BCType type;
  Complex value; // The BC value (if applicable)

  explicit BoundaryCondition(BCType t, Complex v = 0.0) : type(t), value(v) {
  }
};

} // namespace bem::types
