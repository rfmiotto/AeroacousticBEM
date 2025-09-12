#include "domain/physics/GreensFunctions.hpp"
#include "foundation/exceptions/BEMException.hpp"
#include "foundation/utils/Constants.hpp"
#include "foundation/utils/MathUtils.hpp"
#include "types/MathObjectTypes.hpp"
#include <cmath>
#include <stdexcept>

namespace bem::domain::physics {

using bem::foundation::utils::MathUtils::distance;
using bem::foundation::utils::MathUtils::isZero;
using bem::types::Complex;
using bem::types::Point2D;
using bem::types::Real;
using bem::types::Vector2D;
using foundation::utils::Constants::I_4;

Complex hankelFirstKind(int order, const double z) {
  if (isZero(z)) {
    throw bem::foundation::exceptions::BEMException(
        "Hankel function is undefined for z=0.");
  }

  if (order == 0) {
    return {std::cyl_bessel_j(0, z), std::cyl_neumann(0, z)};
  } else if (order == 1) {
    return {std::cyl_bessel_j(1, z), std::cyl_neumann(1, z)};
  } else {
    throw std::invalid_argument(
        "hankelFirstKind is currently only implemented for order 0 and 1.");
  }
}

/**
 * Free-field Green's function G(x,y) for 2D Helmholtz equation
 * G(x,y) = (i/4) * H_0^(1)(k|x-y|)
 * where H_0^(1) is the Hankel function of first kind, order 0
 */
Complex greensFunction2D(const Point2D &source,
                         const Point2D &observer,
                         const double k) {
  if (k <= 0) {
    throw std::invalid_argument("Wave number k must be positive.");
  }

  double r = distance(source, observer);
  if (isZero(r)) {
    throw bem::foundation::exceptions::BEMException(
        "Source and observer points are coincident in greensFunction2D.");
  }
  return I_4 * hankelFirstKind(0, k * r);
}

/**
 * Normal derivative of Green's function
 * ∂G/∂n = (i/4) * k * H_1^(1)(k|x-y|) * (∂r/∂n)
 */
Complex greensFunctionNormalDerivative2D(const Point2D &source,
                                         const Point2D &observer,
                                         const Vector2D &normal,
                                         const double k) {
  if (k <= 0) {
    throw std::invalid_argument("Wave number k must be positive.");
  }

  double r = distance(source, observer);
  if (isZero(r)) {
    throw bem::foundation::exceptions::BEMException(
        "Source and observer points are coincident in "
        "greensFunctionDerivative2D.");
  }

  // Vector from source to observer: diff_vec = observer - source
  Point2D diff_vec = {observer.x - source.x, observer.y - source.y};

  // ∂r/∂n = (r_vec · normal) / r
  // diff_vec / r is the unit vector in the direction of r
  double dr_dn = (diff_vec.x * normal.x + diff_vec.y * normal.y) / r;

  // d/dr H_0^(1)(k*r) = -k * H_1^(1)(k*r)
  // So the full expression is (i/4) * (-k) * H_1^(1)(k*r) * dr/dn
  Complex H1_kr = hankelFirstKind(1, k * r);

  return I_4 * (-k) * H1_kr * dr_dn;
}

} // namespace bem::domain::physics
