#pragma once

#include "foundation/types/BoundaryValueTypes.hpp"

namespace bem::domain::bem {

using types::BCType;
using types::BoundaryCondition;
using types::Real;

/**
 * @brief Interface for boundary condition providers.
 *
 * Implementations provide boundary conditions for a given DOF coordinate.
 */
class IBoundaryConditionProvider {
public:
  virtual ~IBoundaryConditionProvider() = default;

  /// Evaluate BC at a given (x,y) coordinate of a DOF
  [[nodiscard]] virtual BoundaryCondition evaluate(Real x) const = 0;

  /// Type of BC (Dirichlet, Neumann, Kutta, etc.)
  [[nodiscard]] virtual BCType type() const = 0;
};

} // namespace bem::domain::bem
