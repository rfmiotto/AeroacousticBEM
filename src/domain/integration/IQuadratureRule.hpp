#pragma once

#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include <functional>

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Element;
using bem::types::Point2D;
using bem::types::QuadraturePoint;

/**
 * @brief Interface for quadrature rules.
 */
class IQuadratureRule {
public:
  virtual ~IQuadratureRule() = default;

  /**
   * @brief Integrates a function over a panel using the quadrature rule.
   * @param panel The boundary element.
   * @param integrand The function to integrate.
   * @return The result of the integration.
   */
  virtual Complex
  integrate(const Element &panel,
            const std::function<Complex(const Point2D &)> &integrand) const = 0;
};

} // namespace bem::domain::integration
