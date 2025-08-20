#pragma once

#include "domain/integration/IQuadratureRule.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/Constants.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <functional>

namespace bem::domain::integration {

using bem::foundation::utils::Constants::EULER_MASCHERONI;
using bem::foundation::utils::Constants::I_4;
using bem::foundation::utils::Constants::INV_2PI;
using bem::foundation::utils::Constants::PI;

/**
 * @brief Analytical quadrature for constant elements (Kress-like).
 *
 * Implements IQuadratureRule: integrates the single- and double-layer
 * contributions analytically over a constant element.
 *
 * Kress (1987) "On the low wave number asymptotics for the
 * two-dimensional exterior Dirichlet problem for the reduced wave
 * equation"
 */
class AnalyticalConstantQuadrature final : public IQuadratureRule {
public:
  explicit AnalyticalConstantQuadrature(double k) : k_(k) {
  }

  Eigen::ArrayXcd
  integrate(const bem::types::Element &element,
            const std::function<Eigen::ArrayXcd(const QuadraturePoint &,
                                                const bem::types::Element &)>
                & /*unused*/) const override {

    const double xi = 0.0; // center of element in reference [-1,1]
    const Eigen::Index nDOF =
        static_cast<Eigen::Index>(element.shapeFunction(xi).size());

    const double length = element.length();
    if (length <= std::numeric_limits<double>::epsilon()) {
      // Degenerate element: integration domain is zero → return zeros
      single_layer_.setZero(nDOF);
      double_layer_.setZero(nDOF);
      return single_layer_;
    }

    // Single-layer (G)
    const double term =
        (1.0 - EULER_MASCHERONI - std::log(k_ * length * 0.25)) * INV_2PI;
    const std::complex<double> d_val(term * length, 0.25 * length);

    // Double-layer (dG/dn) principal value
    const std::complex<double> h_val(0.5, 0.0);

    single_layer_.setConstant(nDOF, d_val);
    double_layer_.setConstant(nDOF, h_val);

    // By default, return single-layer; assembler can query doubleLayer() if
    // needed
    return single_layer_;
  }

  [[nodiscard]] const bem::types::IntegrationParameters &
  params() const override {
    static bem::types::IntegrationParameters dummy{};
    return dummy;
  }

  const Eigen::ArrayXcd &singleLayer() const {
    return single_layer_;
  }
  const Eigen::ArrayXcd &doubleLayer() const {
    return double_layer_;
  }

private:
  double k_;
  mutable Eigen::ArrayXcd single_layer_;
  mutable Eigen::ArrayXcd double_layer_;
};

} // namespace bem::domain::integration
