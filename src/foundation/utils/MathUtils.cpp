#include "foundation/utils/MathUtils.hpp"
#include "foundation/exceptions/BEMException.hpp" // Required for BEMException
#include "foundation/utils/Constants.hpp" // Required for MACHINE_EPSILON and PI
#include <cmath> // For std::sqrt, std::pow, std::atan2, std::abs, std::cyl_bessel_j, std::cyl_neumann
#include <numeric> // For std::inner_product (potentially, or manual dot product)
#include <stdexcept> // For std::invalid_argument
#include <vector>    // For std::vector

namespace bem::foundation::utils::MathUtils {

double distance(const types::Point2D &p1, const types::Point2D &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return std::sqrt((dx * dx) + (dy * dy));
}

double distanceSquared(const types::Point2D &p1, const types::Point2D &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return (dx * dx) + (dy * dy);
}

double angle(const types::Point2D &p1, const types::Point2D &p2) {
  return std::atan2(p2.y - p1.y, p2.x - p1.x);
}

bool isZero(double value, double tol) { return std::abs(value) < tol; }

bool isEqual(double v1, double v2, double tol) {
  return std::abs(v1 - v2) < tol;
}

double magnitude(types::Complex z) { return std::abs(z); }

double phase(types::Complex z) { return std::arg(z); }

types::Complex hankelFirstKind(int order, double z) {
  if (isZero(z)) {
    // Handle z=0 case carefully, especially for Neumann functions.
    // Y0(0) and Y1(0) are -infinity.
    // Depending on context, might throw, or return specific values if limit is
    // defined. For now, let's rely on std::cyl_neumann behavior or throw if z
    // is zero, as tests use z=1.0, 2.0. The Green's function derivative might
    // call this with k*r. If r is zero, an exception is already thrown there.
    if (order == 0) {
      return {std::cyl_bessel_j(0, z),
              std::cyl_neumann(0, z)}; // Will likely produce inf/nan
    }
    if (order == 1) {
      return {std::cyl_bessel_j(1, z),
              std::cyl_neumann(1, z)}; // Will likely produce inf/nan
    }
    // Fallthrough to throw for other orders or if specific handling for z=0 is
    // needed.
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

types::Complex greensFunction2D(const types::Point2D &source,
                                const types::Point2D &field, double k) {
  double r = distance(source, field);
  if (isZero(r)) {
    throw bem::foundation::exceptions::BEMException(
        "Source and field points are coincident in greensFunction2D.");
  }
  // (i/4) * H_0^(1)(k*r)
  return types::Complex(0.0, 0.25) * hankelFirstKind(0, k * r);
}

types::Complex greensFunctionDerivative2D(const types::Point2D &source,
                                          const types::Point2D &field,
                                          const types::Point2D &normal,
                                          double k) {
  double r = distance(source, field);
  if (isZero(r)) {
    throw bem::foundation::exceptions::BEMException(
        "Source and field points are coincident in "
        "greensFunctionDerivative2D.");
  }

  // Vector from source to field: diff_vec = field - source
  types::Point2D diff_vec = {field.x - source.x, field.y - source.y};

  // dr/dn = dot_product(diff_vec / r, normal)
  // diff_vec / r is the unit vector in the direction of r
  double dr_dn = (diff_vec.x * normal.x + diff_vec.y * normal.y) / r;

  // d/dr H_0^(1)(k*r) = -k * H_1^(1)(k*r)
  // So the full expression is (i/4) * (-k) * H_1^(1)(k*r) * dr/dn
  types::Complex H1_kr = hankelFirstKind(1, k * r);

  return types::Complex(0.0, 0.25) * (-k) * H1_kr * dr_dn;
}

std::pair<std::vector<double>, std::vector<double>>
gaussLegendreQuadrature(int n_points) {
  if (n_points == 2) {
    double point = std::numbers::inv_sqrt3;
    std::vector<double> points = {-point, point};
    std::vector<double> weights = {1.0, 1.0};
    return {points, weights};
  } else {
    throw std::invalid_argument("gaussLegendreQuadrature is currently only "
                                "implemented for n_points = 2.");
  }
}

std::pair<types::Point2D, double>
lineIntegrationTransform(const types::Point2D &start, const types::Point2D &end,
                         double xi) {
  double common_factor_start = (1.0 - xi) / 2.0;
  double common_factor_end = (1.0 + xi) / 2.0;

  types::Point2D transformed_point;
  transformed_point.x =
      common_factor_start * start.x + common_factor_end * end.x;
  transformed_point.y =
      common_factor_start * start.y + common_factor_end * end.y;

  double jacobian = distance(start, end) / 2.0;

  return {transformed_point, jacobian};
}

} // namespace bem::foundation::utils::MathUtils
