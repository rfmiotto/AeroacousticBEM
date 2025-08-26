#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <memory>
#include <omp.h>
#include <vector>

#include "domain/integration/AnalyticalConstantSelfTerm.hpp"
#include "domain/integration/ElementKernelContrib.hpp"
#include "domain/integration/GaussQuadrature.hpp"
#include "domain/integration/IQuadratureRule.hpp"
#include "domain/integration/TellesQuadrature.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/types/MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"

namespace bem::domain::assembly {

using bem::domain::integration::computeVectorizedContribAtQP;
using bem::domain::integration::GaussQuadrature;
using bem::domain::integration::IQuadratureRule;
using bem::domain::integration::KernelType;
using bem::domain::integration::QuadraturePoint;
using bem::domain::integration::TellesQuadrature;

using bem::types::Complex;
using bem::types::Element;
using bem::types::ElementType;
using bem::types::IntegrationParameters;
using bem::types::Point2D;
using bem::types::Real;
using foundation::utils::Constants::GEOMETRY_TOLERANCE;

/**
 * @brief A mapper that, given an element handle (index or ref), returns
 *        the global column indices for that element's local DOFs.
 *
 * Signature:
 *   getGlobalCols(elem_index) -> std::span<const int> (or std::vector<int>)
 *
 * We use a std::function returning std::vector<int> for simplicity.
 */
using GlobalDofMapper = std::function<std::vector<int>(std::size_t elemIndex)>;

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

/**
 * @brief Assemble H and D matrices:
 *        [H]{u} = [D]{q}, for all collocation points.
 *
 * Threading:
 *  - Parallelized over collocation nodes (rows) with OpenMP.
 *  - Quadrature objects are created inside the loop (thread-local) and are
 *    inexpensive (Gauss caches nodes; Telles also caches and just warps).
 *
 * Inputs:
 *  - elements: boundary mesh elements (straight segments).
 *  - collocation: list of collocation points (size = nColloc).
 *  - k: wavenumber.
 *  - ip: integration params (order, distances, flags).
 *  - mapper: maps element index -> global column indices for that element DOFs.
 *
 * Output:
 *  - H, D: complex matrices of size (nColloc x nGlobalDOF).
 */
inline void assembleSystemHD(const std::vector<Element> &elements,
                             const std::vector<Point2D> &collocation,
                             Real k,
                             const IntegrationParameters &ip,
                             const GlobalDofMapper &mapper,
                             Eigen::MatrixXcd &H,
                             Eigen::MatrixXcd &D) {
  const std::size_t nColloc = collocation.size();

  // Determine total number of global DOFs via mapper on last element
  // (Assumes mapper gives valid indices within [0, nGlobalDOF))
  int nGlobal = 0;
  for (std::size_t e = 0; e < elements.size(); ++e) {
    for (int col : mapper(e)) {
      nGlobal = std::max(nGlobal, col + 1);
    }
  }

  H.setZero(static_cast<Eigen::Index>(nColloc), nGlobal);
  D.setZero(static_cast<Eigen::Index>(nColloc), nGlobal);

  // Pre-create a quadrature rule instance usable for far cases (thread-local
  // copy OK) (We still create inside the loop to avoid sharing mutable state
  // across threads.) OpenMP Parallelization over collocation points
  // #pragma omp parallel for schedule(dynamic)
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(nColloc); ++i) {
    const Point2D &x_colloc = collocation[static_cast<std::size_t>(i)];

    // Row views for scatter (thread-safe: different i => different row)
    auto H_row = H.row(i);
    auto D_row = D.row(i);

    for (std::size_t ej = 0; ej < elements.size(); ++ej) {
      const Element &elem = elements[ej];

      // Compute once the column indices for this element
      const std::vector<int> cols = mapper(ej);
      if (cols.empty()) {
        continue;
      }

      // Create small scatter views (RowVectorXcd) into the H/D row
      // covering exactly this element's DOFs.
      Eigen::RowVectorXcd H_scatter(cols.size());
      Eigen::RowVectorXcd D_scatter(cols.size());
      H_scatter.setZero();
      D_scatter.setZero();

      // Select quadrature rule (analytical, Gauss, Telles)
      std::unique_ptr<IQuadratureRule> rule =
          selectQuadratureRule(ip, elem, x_colloc, k);

      // Single-layer (D)
      std::function<Eigen::ArrayXcd(const QuadraturePoint &, const Element &)>
          d_integrand = [&x_colloc, k](const QuadraturePoint &qp,
                                       const Element &e) -> Eigen::ArrayXcd {
        return computeVectorizedContribAtQP(qp, e, x_colloc, k,
                                            KernelType::SINGLE_LAYER_G);
      };
      Eigen::ArrayXcd d_vec = rule->integrate(elem, d_integrand);

      // Double-layer (H)
      std::function<Eigen::ArrayXcd(const QuadraturePoint &, const Element &)>
          h_integrand = [&x_colloc, k](const QuadraturePoint &qp,
                                       const Element &e) -> Eigen::ArrayXcd {
        return computeVectorizedContribAtQP(qp, e, x_colloc, k,
                                            KernelType::DOUBLE_LAYER_DGDN);
      };
      Eigen::ArrayXcd h_vec = rule->integrate(elem, h_integrand);

      // Sanity: match local DOF count
      const auto m = static_cast<Eigen::Index>(cols.size());
      const auto mm = std::min<Eigen::Index>({m, d_vec.size(), h_vec.size()});

      // Use temporary ArrayXcd for safe accumulation
      Eigen::ArrayXcd H_scatter_array =
          Eigen::ArrayXcd::Zero(static_cast<Eigen::Index>(cols.size()));
      Eigen::ArrayXcd D_scatter_array =
          Eigen::ArrayXcd::Zero(static_cast<Eigen::Index>(cols.size()));

      H_scatter_array.head(mm) += h_vec.head(mm);
      D_scatter_array.head(mm) += d_vec.head(mm);

      // Scatter into global row
      for (Eigen::Index c = 0; c < static_cast<Eigen::Index>(cols.size());
           ++c) {
        H_row(cols[c]) += H_scatter_array(c);
        D_row(cols[c]) += D_scatter_array(c);
      }
    }
  }
}

} // namespace bem::domain::assembly
