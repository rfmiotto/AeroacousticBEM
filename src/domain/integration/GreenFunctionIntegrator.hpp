#pragma once

#include "domain/integration/GaussQuadrature.hpp"
#include "domain/integration/SingularityTreatment.hpp"
#include "domain/physics/GreensFunctions.hpp"
#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/Constants.hpp"
#include <memory>

namespace bem::domain::integration {

using bem::domain::physics::greensFunction2D;
using bem::domain::physics::greensFunctionNormalDerivative2D;
using bem::types::Complex;
using bem::types::IntegrationParameters;
using bem::types::Panel;
using bem::types::Point2D;
using bem::types::Vector2D;

class IntegrationStrategy {
public:
  virtual ~IntegrationStrategy() = default;

  [[nodiscard]] virtual Complex
  integrateGreen(const Panel &panel,
                 const Point2D &x,
                 double k,
                 const IntegrationParameters &params) const = 0;

  [[nodiscard]] virtual Complex
  integrateNormalDerivative(const Panel &panel,
                            const Point2D &x,
                            const Vector2D &normal,
                            double k,
                            const IntegrationParameters &params) const = 0;
};

class StandardIntegration : public IntegrationStrategy {
public:
  explicit StandardIntegration(
      int order = bem::foundation::utils::Constants::DEFAULT_INTEGRATION_ORDER,
      std::shared_ptr<SingularityTreatment> singularity_handler =
          std::make_shared<ConstantSingularityTreatment>())
      : quadrature_(order),
        singularity_handler_(std::move(singularity_handler)) {
  }

  [[nodiscard]] Complex
  integrateGreen(const Panel &panel,
                 const Point2D &x,
                 double k,
                 const IntegrationParameters &params) const override {
    if (params.use_singularity_treatment && panel.contains(x)) {
      return singularity_handler_->treatGreen(panel, x, k);
    }

    return quadrature_.integrate(
        panel, [&](const Point2D &y) { return greensFunction2D(x, y, k); });
  }

  [[nodiscard]] Complex integrateNormalDerivative(
      const Panel &panel,
      const Point2D &x,
      const Vector2D &normal,
      double k,
      const IntegrationParameters &params) const override {
    if (params.use_singularity_treatment && panel.contains(x)) {
      return singularity_handler_->treatNormalDerivative(panel, x, normal, k);
    }

    return quadrature_.integrate(panel, [&](const Point2D &y) {
      return greensFunctionNormalDerivative2D(x, y, normal, k);
    });
  }

private:
  GaussQuadrature quadrature_;
  std::shared_ptr<SingularityTreatment> singularity_handler_;
};

} // namespace bem::domain::integration
