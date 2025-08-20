#pragma once

#include <Eigen/Dense>
#include <functional>

#include "domain/integration/GaussLegendreTable.hpp" // QuadraturePoint
#include "foundation/types/GeometryTypes.hpp"        // bem::types::Element
#include "foundation/types/MathObjectTypes.hpp"      // Real, Complex, etc.
#include "types/BEMTypes.hpp"

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Element;
using bem::types::IntegrationParameters;

/**
 * @brief Interface for 1D boundary-element quadrature rules (line segments).
 *
 * The integrand is vectorized: at each quadrature point, it returns an
 * Eigen::ArrayXcd of contributions (e.g., one entry per local shape function).
 * The quadrature rule multiplies by weight and Jacobian and sums, returning a
 * vectorized accumulated result.
 *
 * Design notes:
 * - Geometry (ξ ∈ [-1, 1] → x ∈ ℝ²) is handled by the quadrature via Element.
 * - Jacobian factor is applied inside the quadrature rule.
 * - Implementations (Gauss, Telles, Analytical, etc.) differ only in how
 *   they generate/warp quadrature points/weights.
 */
class IQuadratureRule {
public:
  virtual ~IQuadratureRule() = default;

  /**
   * @brief Vectorized integration over an element.
   *
   * @param element The boundary element to integrate over.
   * @param integrand A function evaluated at each quadrature point that must
   *                  return an ArrayXcd of contributions to be accumulated.
   *                  Signature:
   *                    integrand(const QuadraturePoint& qp,
   *                              const Element& element) -> ArrayXcd
   * @return ArrayXcd Accumulated integral over the element (vectorized).
   */
  virtual Eigen::ArrayXcd
  integrate(const Element &element,
            const std::function<Eigen::ArrayXcd(const QuadraturePoint &,
                                                const Element &)> &integrand)
      const = 0;

  /// @return The integration parameters that configured this rule.
  [[nodiscard]] virtual const IntegrationParameters &params() const = 0;
};

} // namespace bem::domain::integration
