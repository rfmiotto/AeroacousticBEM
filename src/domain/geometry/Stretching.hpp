#pragma once

#include "utils/Constants.hpp"
#include "utils/MathUtils.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace bem::domain::geometry {

using bem::foundation::utils::Constants::GEOMETRY_TOLERANCE;
using bem::foundation::utils::MathUtils::isEqual;

/**
 * @brief Parameters used in the stretching function.
 *
 * These control how the distribution of points is concentrated along the
 * domain.
 *
 * - `paramP`: Controls the linear-to-nonlinear blend.
 * - `paramQ`: Controls the sharpness of the nonlinear clustering.
 * - `concentrateAtEnd`: Controls whether the points must be concentrated
 *                       at the end or beginning of the domain.
 */
struct StretchingParams {
  double paramP = 1.9;
  double paramQ = 2.0;
  bool concentrateAtEnd = true;
};

/**
 * @class StretchingFunction
 * @brief Generates 1D non-uniform node distributions using a smooth stretching
 * function.
 *
 * This class implements a 1D node generation strategy based on a hyperbolic
 * tangent stretching function. It is useful in boundary layer mesh generation
 * and other applications requiring non-uniform point distribution. Here, it
 * is used to properly capture the discontinuities at the leading and trailing
 * edges of the airfoil.
 *
 * The stretching formulation is inspired by:
 * - G.O. Roberts, "Computational meshes for boundary layer problems", in: M.
 * Holt (Ed.), Numerical Methods in Fluid Dynamics, Lecture Notes in Physics,
 * vol. 8, Springer, 1971, pp. 171–177.
 * - P.R. Eiseman, "A multi-surface method of coordinate generation", J. Comput.
 * Phys. 33 (1979) 118–150.
 */
class StretchingFunction {
public:
  /**
   * @brief Constructs a StretchingFunction over a given domain with optional
   * parameters.
   *
   * @param xBeg Starting coordinate of the domain.meaning
   * @param xEnd Ending coordinate of the domain.
   * @param numPoints Number of nodes to generate (must be >= 2).
   * @param params Stretching parameters (optional; default values provided).
   *
   * @throws std::invalid_argument if xEnd <= xBeg or numPoints < 2.
   */
  StretchingFunction(double xBeg,
                     double xEnd,
                     int numPoints,
                     const StretchingParams &params = StretchingParams{})
      : paramP_(params.paramP), paramQ_(params.paramQ), xBeg_(xBeg),
        xEnd_(xEnd), numPoints_(numPoints),
        concentrateAtEnd_(params.concentrateAtEnd) {

    if (numPoints_ < 2) {
      throw std::invalid_argument("numPoints must be at least 2.");
    }

    if (isEqual(xBeg, xEnd, GEOMETRY_TOLERANCE)) {
      throw std::invalid_argument("numPoints must be at least 2.");
    }
  }

  /**
   * @brief Generates the non-uniformly spaced 1D grid.
   *
   * The stretching function clusters points according to the provided
   * parameters and orientation.
   *
   * @return A vector of node coordinates between xBeg and xEnd.
   */
  [[nodiscard]] std::vector<double> generateNodes() const {
    std::vector<double> xnodes;
    xnodes.reserve(numPoints_);

    // Stretching function for parameter s
    std::vector<double> sValues;
    sValues.reserve(numPoints_);

    for (int i = 0; i < numPoints_; ++i) {
      double eta = static_cast<double>(i) / (numPoints_ - 1);
      double s = (paramP_ * eta) +
                 ((1.0 - paramP_) * (1.0 - std::tanh(paramQ_ * (1.0 - eta)) /
                                               std::tanh(paramQ_)));
      sValues.push_back(s);
    }

    // Determine direction and concentration logic
    if (!concentrateAtEnd_) {
      for (auto &s : sValues) {
        s = 1.0 - s;
      }
      std::ranges::reverse(sValues);
    }

    // Compute the physical coordinates
    for (double s : sValues) {
      double x = xBeg_ + (s * (xEnd_ - xBeg_));
      xnodes.push_back(x);
    }
    return xnodes;
  }

  /**
   * @brief Generates the analytical derivatives of the node coordinates with
   * respect to the normalized coordinate eta.
   *
   * @return A vector of derivatives dx/deta at each node.
   */
  [[nodiscard]] std::vector<double> generateNodeDerivatives() const {
    std::vector<double> derivatives;
    derivatives.reserve(numPoints_);

    const double tanhQ = std::tanh(paramQ_);

    for (int i = 0; i < numPoints_; ++i) {
      double eta = static_cast<double>(i) / (numPoints_ - 1);
      double sechSquared = 1.0 - std::pow(std::tanh(paramQ_ * (1.0 - eta)), 2);
      double ds_deta =
          paramP_ + ((1.0 - paramP_) * paramQ_ * sechSquared / tanhQ);
      derivatives.push_back(ds_deta);
    }

    // Determine direction and concentration logic
    if (!concentrateAtEnd_) {
      for (auto &ds_deta : derivatives) {
        ds_deta *= -1.0;
      }
      std::ranges::reverse(derivatives);
    }

    // Compute the physical derivatives: dx/deta = ds/deta * (xEnd - xBeg)
    const double dx = xEnd_ - xBeg_;
    for (auto &ds_deta : derivatives) {
      ds_deta *= dx;
    }

    return derivatives;
  }

private:
  double paramP_;         ///< Blending parameter (linear ↔ nonlinear)
  double paramQ_;         ///< Stretching sharpness
  double xBeg_;           ///< Start of domain
  double xEnd_;           ///< End of domain
  int numPoints_;         ///< Number of points to generate
  bool concentrateAtEnd_; ///< Extreme where to concentrate the points
};

} // namespace bem::domain::geometry
