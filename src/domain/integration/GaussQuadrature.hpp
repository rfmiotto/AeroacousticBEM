#pragma once

#include <Eigen/Dense>
#include <functional>

#include "domain/integration/GaussLegendreTable.hpp"
#include "domain/integration/IQuadratureRule.hpp"
#include "foundation/types/GeometryTypes.hpp"

namespace bem::domain::integration {

/**
 * @brief Gauss-Legendre quadrature rule for boundary element integration.
 *
 * This class caches the reference quadrature nodes (ξ, w) once at construction
 * and maps them to physical coordinates on each integrate() call.
 */
class GaussQuadrature final : public IQuadratureRule {
public:
  explicit GaussQuadrature(const bem::types::IntegrationParameters &ip)
      : params_(ip), qps_ref_(GaussLegendreTable::getPoints(ip.order)) {
  }

  Eigen::ArrayXcd
  integrate(const bem::types::Element &element,
            const std::function<Eigen::ArrayXcd(const QuadraturePoint &,
                                                const bem::types::Element &)>
                &integrand) const override {
    using bem::types::Real;

    // Jacobian for the mapping from ξ ∈ [-1,1] to element (constant for
    // straight elements)
    const Real J = element.jacobian();

    Eigen::ArrayXcd acc; // lazy allocation
    bool initialized = false;

    for (const auto &qp_ref : qps_ref_) {
      // Map reference ξ to physical point on the element
      const auto x = element.parametricPoint(qp_ref.xi);

      // Quadrature point with updated physical location
      QuadraturePoint qp{qp_ref.xi, qp_ref.weight, x};

      // Evaluate integrand (vectorized contributions, e.g. per shape fn)
      Eigen::ArrayXcd contrib = integrand(qp, element);

      if (!initialized) {
        acc = Eigen::ArrayXcd::Zero(contrib.size());
        initialized = true;
      }

      // Accumulate: weight * Jacobian * contribution
      acc += (qp.weight * J) * contrib;
    }

    if (!initialized) {
      // Edge case: zero order? (should not happen) Return empty.
      return Eigen::ArrayXcd{}; // safety: no quadrature points
    }
    return acc;
  }

  [[nodiscard]] const bem::types::IntegrationParameters &
  params() const override {
    return params_;
  }

private:
  bem::types::IntegrationParameters params_;

  // Cached Gauss-Legendre reference quadrature nodes (ξ, w).
  // Physical mapping is applied per integrate() call.
  std::vector<QuadraturePoint> qps_ref_;
};

} // namespace bem::domain::integration
