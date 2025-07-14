#pragma once

#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/Constants.hpp"
#include <cmath>

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Element;
using bem::types::Point2D;
using bem::types::Vector2D;

/* Kress (1987) "On the low wave number asymptotics for the
 * two-dimensional exterior Dirichlet problem for the
 * reduced wave equation"
 */
class SingularityTreatment {
public:
  virtual ~SingularityTreatment() = default;

  [[nodiscard]] virtual Complex
  treatGreen(const Element &element, const Point2D &x, double k) const = 0;

  [[nodiscard]] virtual Complex treatNormalDerivative(const Element &element,
                                                      const Point2D &x,
                                                      const Vector2D &normal,
                                                      double k) const = 0;
};

using bem::foundation::utils::Constants::EULER_MASCHERONI;
using bem::foundation::utils::Constants::I_4;
using bem::foundation::utils::Constants::INV_2PI;
using bem::foundation::utils::Constants::PI;

class ConstantSingularityTreatment : public SingularityTreatment {
public:
  [[nodiscard]] Complex treatGreen(const Element &element,
                                   const Point2D & /*x*/,
                                   double k) const override {
    const double length = element.length();
    const double term =
        (1.0 - EULER_MASCHERONI - std::log(k * length * 0.25)) * INV_2PI;
    return Complex(term, 0.25) * length;
  }

  [[nodiscard]] Complex treatNormalDerivative(const Element & /*panel*/,
                                              const Point2D & /*x*/,
                                              const Vector2D & /*normal*/,
                                              double /*k*/) const override {
    return 0.5;
  }
};

} // namespace bem::domain::integration
