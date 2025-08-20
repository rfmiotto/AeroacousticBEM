#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "domain/integration/GaussLegendreTable.hpp" // QuadraturePoint
#include "domain/integration/IQuadratureRule.hpp"    // Interface
#include "foundation/types/GeometryTypes.hpp" // Element, Point2D, Vector2D

namespace bem::domain::integration {

using bem::types::Element;
using bem::types::IntegrationParameters;
using bem::types::Point2D;
using bem::types::Real;
using bem::types::Vector2D;

/**
 * @brief Telles quadrature rule for near-singular integrals on straight
 *        boundary elements.
 *
 * Pattern mirrors GaussQuadrature:
 *  - Caches reference Gauss nodes/weights (s_k, w_k) in the ctor.
 *  - On integrate(), computes the singular-attractor ξ* from collocation x*,
 *    warps (s_k, w_k) -> (ξ_k, W_k) using a Telles-like transformation,
 *    maps ξ_k -> x_k, and accumulates vectorized integrand with element
 * Jacobian.
 *
 * IMPORTANT: The actual Telles mapping is encapsulated in tellesWarp(). You can
 * replace that body with your preferred/validated Telles formulation while
 * keeping the rest intact.
 */
class TellesQuadrature final : public IQuadratureRule {
public:
  /**
   * @param ip         Integration parameters (order, tolerances, etc.).
   * @param x_star     Collocation / observer point that induces
   * near-singularity. The node clustering will be focused toward the projection
   *                   of x_star onto the current element (computed per
   * integrate()).
   */
  TellesQuadrature(const IntegrationParameters &ip, const Point2D &x_star)
      : params_(ip), x_star_(x_star),
        s_w_ref_(GaussLegendreTable::getPoints(ip.order)) {
  }

  [[nodiscard]] const IntegrationParameters &params() const override {
    return params_;
  }

  /**
   * @brief Vectorized integration using Telles warping.
   *
   * The integrand must return ArrayXcd (e.g., one entry per local shape fn).
   * We compute ξ* from the collocation x_star_ and the current element, warp
   * the cached (s, w), then accumulate: sum_k [ W_k * J_elem * integrand(qp_k)
   * ].
   */
  Eigen::ArrayXcd
  integrate(const Element &element,
            const std::function<Eigen::ArrayXcd(const QuadraturePoint &,
                                                const Element &)> &integrand)
      const override {
    const Real J_elem = element.jacobian(); // length/2 for straight segments

    // 1) Compute ξ* (projection of x_star_ onto element, mapped to [-1,1])
    const Real xi_star = projectToXi(element, x_star_);

    // 2) Prepare accumulator
    Eigen::ArrayXcd acc;
    bool init = false;

    // 3) Loop over cached reference Gauss nodes/weights (s, w)
    for (const auto &qp_ref : s_w_ref_) {
      // 3a) Telles warp: s -> ξ, with warped weight W = w * dξ/ds
      Real xi_warp = 0.0;
      Real w_warp = 0.0;
      tellesWarp(qp_ref.xi, qp_ref.weight, xi_star, xi_warp, w_warp);

      // Safety clamp to [-1,1] (numerical): keep mapping valid
      xi_warp = std::max(xi_warp, -1.0);
      xi_warp = std::min(xi_warp, 1.0);

      // 3b) Map to physical point on the element
      const auto x_phys = element.parametricPoint(xi_warp);

      // Build quadrature point with *warped* (ξ, W) and physical coordinate
      QuadraturePoint qp{xi_warp, w_warp, x_phys};

      // 3c) Evaluate vectorized integrand at this quadrature point
      Eigen::ArrayXcd contrib = integrand(qp, element);

      if (!init) {
        acc = Eigen::ArrayXcd::Zero(contrib.size());
        init = true;
      }

      // 3d) Accumulate with element Jacobian
      acc += (qp.weight * J_elem) * contrib;
    }

    return init ? acc : Eigen::ArrayXcd{};
  }

private:
  IntegrationParameters params_;
  Point2D x_star_; // observer/collocation (near-singularity attractor)

  // Cached reference Gauss nodes/weights: s \in [-1,1], weight w
  std::vector<QuadraturePoint> s_w_ref_;

  // --- Helpers --------------------------------------------------------------

  /**
   * @brief Project physical point x_star onto the straight element and return
   *        the corresponding ξ ∈ [-1,1].
   *
   * For a straight segment with endpoints A (ξ=-1) and B (ξ=+1), we:
   *  - Compute the param t ∈ [0,1] along AB for the orthogonal projection of
   * x_star onto the infinite line through AB, clamped to [0,1] for the segment.
   *  - Map t -> ξ = 2t - 1.
   */
  static Real projectToXi(const Element &elem, const Point2D &x_star) {
    const auto &A = elem.start;
    const auto &B = elem.end;

    const Real ax = A.x;
    const Real ay = A.y;
    const Real bx = B.x;
    const Real by = B.y;

    const Real vx = bx - ax;
    const Real vy = by - ay;
    const Real wx = x_star.x - ax;
    const Real wy = x_star.y - ay;

    const Real vv = (vx * vx) + (vy * vy);
    Real t = (vv > Real(0)) ? (wx * vx + wy * vy) / vv
                            : Real(0.5); // center if degenerate
    t = std::max(t, Real(0));
    t = std::min(t, Real(1));

    // Map t in [0,1] to ξ in [-1,1]
    return (Real(2) * t) - Real(1);
  }

  /**
   * @brief Telles-like cubic warping that clusters nodes toward ξ_star.
   *
   * INPUT:
   *   s     : reference Gauss node in [-1,1]
   *   w     : reference Gauss weight (>0)
   *   xiStar: target in [-1,1] that attracts nodes (near-singular location)
   * OUTPUT:
   *   xi    : warped abscissa in [-1,1]
   *   W     : warped weight (= w * dξ/ds at s)
   *
   * Implementation note:
   *   This function is intentionally isolated. If you have your canonical
   *   Telles transform (per your preferred reference), drop it here and keep
   *   the rest of the class unchanged. The current implementation is a smooth
   *   cubic clustering that:
   *     - preserves endpoints (s=±1 -> ξ=±1),
   *     - preserves symmetry if xiStar=0,
   *     - concentrates nodes near xiStar,
   *     - provides an analytic derivative for weight update.
   *
   *   Mapping used:
   *     Let a = xiStar. Define blend parameter p∈[0,1) based on proximity of a
   * to the element (closer -> stronger clustering). We use: p = 0.85 * (0.5 +
   * 0.5 * |a|)   // heuristic, monotone in |a| Define a cubic that preserves
   * endpoints and slopes at s=±1 approximately: ξ(s) = s + p * ( (1 - s*s) * (a
   * - s) ) Then: dξ/ds = 1 + p * ( -2s*(a - s) - (1 - s*s) )
   *
   *   You can replace ξ(s) and dξ/ds below with your validated Telles mapping.
   */
  static void tellesWarp(Real s, Real w, Real xiStar, Real &xi, Real &W) {
    // Strength of clustering (heuristic): stronger when |xiStar| near 1 (closer
    // to element ends)
    const Real p = Real(0.85) * (Real(0.5) + Real(0.5) * std::abs(xiStar));

    // Smooth cubic clustering toward xiStar; preserves endpoints; analytic
    // derivative
    const Real one_minus_ss = Real(1) - (s * s);
    const Real xi_map = s + (p * (one_minus_ss * (xiStar - s)));

    // Derivative dξ/ds for weight update
    const Real dxi_ds =
        Real(1) + (p * (-Real(2) * s * (xiStar - s) - one_minus_ss));

    xi = xi_map;
    W = w * std::abs(dxi_ds);
  }
};

} // namespace bem::domain::integration
