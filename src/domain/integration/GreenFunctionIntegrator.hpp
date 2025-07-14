#pragma once

#include "domain/integration/GaussQuadrature.hpp"
#include "domain/integration/SingularityTreatment.hpp"
#include "domain/integration/TellesQuadrature.hpp"
#include "domain/physics/GreensFunctions.hpp"
#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"

#include <functional>
#include <memory>

namespace bem::domain::integration {

using bem::domain::physics::greensFunction2D;
using bem::domain::physics::greensFunctionNormalDerivative2D;
using bem::types::Complex;
using bem::types::Element;
using bem::types::ElementType;
using bem::types::IntegrationParameters;
using bem::types::Point2D;
using bem::types::Vector2D;

using GreenIntegrator =
    std::function<Complex(const Element &, const Point2D &, double)>;
using NormalDerivativeIntegrator = std::function<
    Complex(const Element &, const Point2D &, const Vector2D &, double)>;

inline GreenIntegrator makeGreenIntegrator(
    const IntegrationParameters &params,
    const std::shared_ptr<SingularityTreatment> &singularity_handler =
        std::make_shared<ConstantSingularityTreatment>()) {
  const auto order = params.order;

  return [=](const Element &element, const Point2D &x, double k) -> Complex {
    const bool is_singular = element.contains(x);

    if (is_singular && element.type == ElementType::CONSTANT) {
      return singularity_handler->treatGreen(element, x, k);
    }

    if (is_singular && params.use_singularity_treatment) {
      TellesQuadrature telles(order);
      return telles.integrate(element, x, [&](const Point2D &y) {
        return greensFunction2D(x, y, k);
      });
    }

    // Regular integration
    GaussQuadrature gauss(order);
    return gauss.integrate(element, x, [&](const Point2D &y) {
      return greensFunction2D(x, y, k);
    });
  };
}

inline NormalDerivativeIntegrator makeNormalDerivativeIntegrator(
    const IntegrationParameters &params,
    const std::shared_ptr<SingularityTreatment> &singularity_handler =
        std::make_shared<ConstantSingularityTreatment>()) {
  const auto order = params.order;

  return [=](const Element &element, const Point2D &x, const Vector2D &normal,
             double k) -> Complex {
    const bool is_singular = element.contains(x);

    if (is_singular && element.type == ElementType::CONSTANT) {
      return singularity_handler->treatNormalDerivative(element, x, normal, k);
    }

    if (is_singular && params.use_singularity_treatment) {
      TellesQuadrature telles(order);
      return telles.integrate(element, x, [&](const Point2D &y) {
        return greensFunctionNormalDerivative2D(x, y, normal, k);
      });
    }

    GaussQuadrature gauss(order);
    return gauss.integrate(element, x, [&](const Point2D &y) {
      return greensFunctionNormalDerivative2D(x, y, normal, k);
    });
  };
}

} // namespace bem::domain::integration
