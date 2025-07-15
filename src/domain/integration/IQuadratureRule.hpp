#pragma once

#include "GaussLegendreTable.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include <functional>

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Element;
using bem::types::Point2D;

/**
 * @brief Interface for quadrature rules.
 */
class IQuadratureRule {
public:
  virtual ~IQuadratureRule() = default;

  /**
   * @brief Integrates a function over an element using the quadrature rule.
   * @param element The boundary element.
   * @param integrand The function to integrate.
   * @return The result of the integration.
   */
  virtual Complex
  integrate(const Element &element,
            const Point2D &x,
            const std::function<Complex(const Point2D &)> &integrand) const = 0;
};

} // namespace bem::domain::integration
