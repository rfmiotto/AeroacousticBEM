#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <memory>
#include <omp.h>

#include "domain/integration/AnalyticalConstantSelfTerm.hpp"
#include "domain/integration/ElementKernelContrib.hpp"
#include "domain/integration/GaussQuadrature.hpp"
#include "domain/integration/IQuadratureRule.hpp"
#include "domain/integration/TellesQuadrature.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/types/MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"

namespace bem::domain::integration {

using bem::domain::integration::computeVectorizedContribAtQP;
using bem::domain::integration::GaussQuadrature;
using bem::domain::integration::IQuadratureRule;
using bem::domain::integration::TellesQuadrature;

using bem::types::Element;
using bem::types::ElementType;
using bem::types::IntegrationParameters;
using bem::types::Point2D;
using bem::types::Real;
using foundation::utils::Constants::GEOMETRY_TOLERANCE;

/**
 * @brief Simple near-singularity predicate: distance from collocation to
 *        element segment < params.singularity_distance.
 *        You can replace with your more sophisticated test if desired.
 */
inline bool
isNearSingular(const Element &elem, const Point2D &x, Real threshold) {
  // Distance from point to segment (2D)
  const auto &A = elem.start;
  const auto &B = elem.end;

  const Real vx = B.x - A.x;
  const Real vy = B.y - A.y;
  const Real wx = x.x - A.x;
  const Real wy = x.y - A.y;

  const Real vv = (vx * vx) + (vy * vy);
  Real t = (vv > Real(0)) ? (wx * vx + wy * vy) / vv : Real(0.5);
  t = std::max(t, Real(0));
  t = std::min(t, Real(1));

  const Real px = A.x + (t * vx);
  const Real py = A.y + (t * vy);

  const Real dx = x.x - px;
  const Real dy = x.y - py;
  const Real dist = std::sqrt((dx * dx) + (dy * dy));
  return dist <= threshold;
}

/**
 * @brief Strategy for selecting which quadrature to use for a given
 *        collocation point and element.
 *
 *  - If collocation lies on constant element -> Analytical self term.
 *  - Else if near singular -> Telles
 *  - Else -> Gauss
 *
 * Returns:
 *   - std::nullopt for Analytical (handled separately),
 *   - otherwise a unique_ptr<IQuadratureRule> to use.
 */
inline std::unique_ptr<IQuadratureRule>
selectQuadratureRule(const IntegrationParameters &ip,
                     const Element &elem,
                     const Point2D &x_colloc,
                     double k) {
  // Collocation lies on element & element is constant -> analytical
  if (elem.type == ElementType::CONSTANT &&
      elem.contains(x_colloc, GEOMETRY_TOLERANCE)) {
    return std::make_unique<integration::AnalyticalConstantQuadrature>(k);
  }

  // Near-singular -> Telles
  if (isNearSingular(elem, x_colloc, ip.singularity_distance)) {
    return std::make_unique<TellesQuadrature>(ip, x_colloc);
  }

  // Default -> Gauss
  return std::make_unique<GaussQuadrature>(ip);
}

} // namespace bem::domain::integration
