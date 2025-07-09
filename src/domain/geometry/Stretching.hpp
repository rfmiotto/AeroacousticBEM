#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * @brief Parameters used in the stretching function.
 *
 * These control how the distribution of points is concentrated along the
 * domain.
 *
 * - `paramP`: Controls the linear-to-nonlinear blend.
 * - `paramQ`: Controls the sharpness of the nonlinear clustering.
 */
struct StretchingParams {
  double paramP = 1.9;
  double paramQ = 2.0;
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
   * @param xBeg Starting coordinate of the domain.
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
        xEnd_(xEnd), numPoints_(numPoints) {

    if (xEnd_ <= xBeg_) {
      throw std::invalid_argument("xEnd must be greater than xBeg.");
    }
    if (numPoints_ < 2) {
      throw std::invalid_argument("numPoints must be at least 2.");
    }
  }

  /**
   * @brief Generates the non-uniformly spaced 1D grid.
   *
   * The stretching function clusters points according to the provided
   * parameters.
   *
   * @return A vector of node coordinates between xBeg and xEnd.
   */
  [[nodiscard]] std::vector<double> generateNodes() const {
    std::vector<double> xnodes;
    xnodes.reserve(numPoints_);

    for (int i = 0; i < numPoints_; ++i) {
      double eta = static_cast<double>(i) / (numPoints_ - 1);
      double s = (paramP_ * eta) +
                 ((1.0 - paramP_) * (1.0 - std::tanh(paramQ_ * (1.0 - eta)) /
                                               std::tanh(paramQ_)));
      double x = xBeg_ + ((xEnd_ - xBeg_) * s);
      xnodes.push_back(x);
    }

    return xnodes;
  }

private:
  double paramP_; ///< Blending parameter (linear ↔ nonlinear)
  double paramQ_; ///< Stretching sharpness
  double xBeg_;   ///< Start of domain
  double xEnd_;   ///< End of domain
  int numPoints_; ///< Number of points to generate
};
