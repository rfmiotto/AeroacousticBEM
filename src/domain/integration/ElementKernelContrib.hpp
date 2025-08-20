#pragma once

#include <Eigen/Dense>
#include <vector>

#include "domain/integration/GaussLegendreTable.hpp"
#include "domain/physics/GreensFunctions.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/types/MathObjectTypes.hpp"

namespace bem::domain::integration {

using bem::types::Complex;
using bem::types::Element;
using bem::types::Point2D;
using bem::types::Real;

/**
 * @brief Compute vectorized contributions at one quadrature point for either:
 *        - Single-layer kernel (G):   contributes to D-matrix (relates flux q)
 *        - Double-layer kernel (∂G/∂n): contributes to H-matrix (relates u)
 *
 * Returns an ArrayXcd with size = number of local shape functions on `element`.
 * Each entry corresponds to the contribution for that local basis function,
 * ready to be scattered into the global matrix columns.
 *
 * This function is meant to be used as the `integrand` callback for
 * IQuadratureRule::integrate(…).
 */

// Which kernel to evaluate at the quadrature point
enum class KernelType : std::uint8_t { SINGLE_LAYER_G, DOUBLE_LAYER_DGDN };

inline Eigen::ArrayXcd computeVectorizedContribAtQP(const QuadraturePoint &qp,
                                                    const Element &element,
                                                    const Point2D &collocation,
                                                    Real wavenumber_k,
                                                    KernelType kernel) {
  // Evaluate local shape functions at ξ
  const std::vector<Real> N = element.shapeFunction(qp.xi);

  // Map std::vector -> Eigen::ArrayXd, then cast to Complex for scaling
  Eigen::Map<const Eigen::ArrayXd> N_arr(N.data(),
                                         static_cast<Eigen::Index>(N.size()));
  const Eigen::ArrayXcd N_c =
      N_arr.cast<Complex>(); // broadcast-friendly complex view

  // Kernels:
  //   G(x,y) with source at boundary point qp.point, observer at collocation
  //   ∂G/∂n at qp.point with element's outward normal
  Complex scalar{};
  if (kernel == KernelType::SINGLE_LAYER_G) {
    scalar = bem::domain::physics::greensFunction2D(qp.point, collocation,
                                                    wavenumber_k);
  } else {
    const auto n = element.normal(); // outward normal at the (straight) element
    scalar = bem::domain::physics::greensFunctionNormalDerivative2D(
        qp.point, collocation, n, wavenumber_k);
  }

  // Vectorized contribution: scalar kernel * shape functions (per local DOF)
  // (Jacobian and quadrature weight are applied by the quadrature rule.)
  return (scalar * N_c);
}

} // namespace bem::domain::integration
