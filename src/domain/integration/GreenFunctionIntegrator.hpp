#pragma once

#include "domain/integration/GaussQuadrature.hpp"
#include "domain/physics/GreensFunctions.hpp"
#include "foundation/exceptions/BEMException.hpp"
#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/Constants.hpp"
#include "foundation/utils/MathUtils.hpp"

namespace bem::domain::integration {

using bem::domain::physics::greensFunction2D;
using bem::domain::physics::greensFunctionNormalDerivative2D;
using bem::foundation::exceptions::BEMException;
using bem::foundation::utils::Constants::I_4;
using bem::foundation::utils::Constants::PI;
using bem::foundation::utils::Constants::SINGULARITY_DISTANCE_THRESHOLD;
using bem::foundation::utils::MathUtils::distance;
using bem::types::Complex;
using bem::types::IntegrationParameters;
using bem::types::LineSegment;
using bem::types::Point2D;
using bem::types::Vector2D;

class IntegrationStrategy {
public:
  virtual ~IntegrationStrategy() = default;

  [[nodiscard]] virtual Complex
  integrateGreen(const LineSegment &panel,
                 const Point2D &x,
                 double k,
                 const IntegrationParameters &params) const = 0;

  [[nodiscard]] virtual Complex
  integrateNormalDerivative(const LineSegment &panel,
                            const Point2D &x,
                            const Vector2D &normal,
                            double k,
                            const IntegrationParameters &params) const = 0;
};

class StandardIntegration : public IntegrationStrategy {
public:
  explicit StandardIntegration(
      int order = bem::foundation::utils::Constants::DEFAULT_INTEGRATION_ORDER)
      : quadrature_(order) {
  }

  [[nodiscard]] Complex
  integrateGreen(const LineSegment &panel,
                 const Point2D &x,
                 double k,
                 const IntegrationParameters &params) const override {
    double dist = distance(x, panel.midpoint());
    if (params.use_singularity_treatment &&
        dist < params.singularity_distance) {
      double length = panel.length();
      double a = k * length / 2;
      return -I_4 * length / (2 * PI) * (std::log(a) - 1); // FIXME: check
    }

    return quadrature_.integrate(
        panel, [&x, k](const Point2D &y) { return greensFunction2D(x, y, k); });
  }

  [[nodiscard]] Complex integrateNormalDerivative(
      const LineSegment &panel,
      const Point2D &x,
      const Vector2D &normal,
      double k,
      const IntegrationParameters &params) const override {
    double dist = distance(x, panel.midpoint());
    if (params.use_singularity_treatment &&
        dist < params.singularity_distance) {
      return 0.5;
    }

    return quadrature_.integrate(panel, [&x, &normal, k](const Point2D &y) {
      return greensFunctionNormalDerivative2D(x, y, normal, k);
    });
  }

private:
  GaussQuadrature quadrature_;
};

} // namespace bem::domain::integration
