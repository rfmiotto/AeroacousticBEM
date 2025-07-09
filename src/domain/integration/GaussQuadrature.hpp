#pragma once

#include "foundation/types/MathObjectTypes.hpp"
#include <array>
#include <vector>

namespace bem::domain {

// Integration point for numerical quadrature
struct IntegrationPoint {
  types::Point2D coordinates;
  double weight;

  IntegrationPoint(const types::Point2D &coords, double w)
      : coordinates(coords), weight(w) {
  }
};

/**
 * Gauss-Legendre quadrature rules for numerical integration
 */
class GaussQuadrature {
public:
  explicit GaussQuadrature(int order) : order_(order) {
    setupQuadraturePoints();
  }

  /**
   * Get integration points and weights for interval [-1, 1]
   */
  const std::vector<IntegrationPoint> &getPoints() const {
    return points_;
  }

  /**
   * Transform integration points from [-1, 1] to [a, b]
   */
  std::vector<IntegrationPoint> transformInterval(double a, double b) const {
    std::vector<IntegrationPoint> transformed;
    transformed.reserve(points_.size());

    double jacobian = 0.5 * (b - a);
    double midpoint = 0.5 * (a + b);

    for (const auto &point : points_) {
      double xi = jacobian * point.coordinates.x + midpoint;
      double weight = jacobian * point.weight;
      transformed.emplace_back(types::Point2D(xi, 0.0), weight);
    }

    return transformed;
  }

  /**
   * Integrate over a line element from point1 to point2
   */
  template <typename Func>
  types::Complex integrateLineElement(const types::Point2D &point1,
                                      const types::Point2D &point2,
                                      Func integrand) const {
    types::Vector2D tangent = point2 - point1;
    double length = tangent.norm();

    if (length < constants::EPS) {
      return Complex(0.0, 0.0);
    }

    tangent = tangent * (1.0 / length); // normalize

    types::Complex result(0.0, 0.0);

    // Transform to [0, 1] parameter space
    for (const auto &point : points_) {
      double xi = 0.5 * (point.coordinates.x + 1.0); // [-1,1] -> [0,1]
      types::Point2D eval_point = point1 + tangent * (xi * length);

      types::Complex value = integrand(eval_point);
      result += value * point.weight;
    }

    return result * (0.5 * length); // Jacobian for [0,1] transformation
  }

private:
  int order_;
  std::vector<IntegrationPoint> points_;

  void setupQuadraturePoints() {
    switch (order_) {
    case 1:
      setupGauss1();
      break;
    case 2:
      setupGauss2();
      break;
    case 3:
      setupGauss3();
      break;
    case 4:
      setupGauss4();
      break;
    case 8:
      setupGauss8();
      break;
    default:
      throw BEMException("Unsupported quadrature order: " +
                         std::to_string(order_));
    }
  }

  void setupGauss1() {
    points_.clear();
    points_.emplace_back(types::Point2D(0.0, 0.0), 2.0);
  }

  void setupGauss2() {
    points_.clear();
    double xi = 1.0 / std::sqrt(3.0);
    points_.emplace_back(types::Point2D(-xi, 0.0), 1.0);
    points_.emplace_back(types::Point2D(xi, 0.0), 1.0);
  }

  void setupGauss3() {
    points_.clear();
    double xi = std::sqrt(3.0 / 5.0);
    points_.emplace_back(types::Point2D(-xi, 0.0), 5.0 / 9.0);
    points_.emplace_back(types::Point2D(0.0, 0.0), 8.0 / 9.0);
    points_.emplace_back(types::Point2D(xi, 0.0), 5.0 / 9.0);
  }

  void setupGauss4() {
    points_.clear();
    double xi1 = std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
    double xi2 = std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
    double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
    double w2 = (18.0 - std::sqrt(30.0)) / 36.0;

    points_.emplace_back(types::Point2D(-xi2, 0.0), w2);
    points_.emplace_back(types::Point2D(-xi1, 0.0), w1);
    points_.emplace_back(types::Point2D(xi1, 0.0), w1);
    points_.emplace_back(types::Point2D(xi2, 0.0), w2);
  }

  void setupGauss8() {
    points_.clear();
    // 8-point Gauss-Legendre quadrature
    std::array<double, 4> xi = {0.1834346424956498, 0.5255324099163290,
                                0.7966664774136267, 0.9602898564975363};

    std::array<double, 4> w = {0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};

    for (int i = 0; i < 4; ++i) {
      points_.emplace_back(types::Point2D(-xi[i], 0.0), w[i]);
      points_.emplace_back(types::Point2D(xi[i], 0.0), w[i]);
    }
  }
};

/**
 * Adaptive integration for singular and near-singular integrals
 */
class AdaptiveIntegration {
public:
  AdaptiveIntegration(double tolerance = 1e-6, int max_levels = 10)
      : tolerance_(tolerance), max_levels_(max_levels) {
  }

  /**
   * Adaptive integration over line element with singularity handling
   */
  template <typename Func>
  types::Complex
  integrateAdaptive(const types::Point2D &point1,
                    const types::Point2D &point2,
                    Func integrand,
                    const types::Point2D &singularity_point) const {
    // Check if singularity is close to the integration interval
    double dist_to_line =
        distancePointToLine(singularity_point, point1, point2);

    if (dist_to_line < tolerance_ * (point2 - point1).norm()) {
      return integrateSingular(point1, point2, integrand, singularity_point);
    } else {
      // Use standard Gauss quadrature
      GaussQuadrature gauss(8);
      return gauss.integrateLineElement(point1, point2, integrand);
    }
  }

private:
  double tolerance_;
  int max_levels_;

  double distancePointToLine(const types::Point2D &point,
                             const types::Point2D &line_start,
                             const types::Point2D &line_end) const {
    Vector2D line_vec = line_end - line_start;
    Vector2D point_vec = point - line_start;

    double line_length = line_vec.norm();
    if (line_length < constants::EPS) {
      return (point - line_start).norm();
    }

    double projection = point_vec.dot(line_vec) / line_length;
    projection = std::max(0.0, std::min(line_length, projection));

    Point2D closest_point = line_start + line_vec * (projection / line_length);
    return (point - closest_point).norm();
  }

  template <typename Func>
  types::Complex
  integrateSingular(const types::Point2D &point1,
                    const types::Point2D &point2,
                    Func integrand,
                    const types::Point2D &singularity_point) const {
    // Implementation of singular integration techniques
    // This is a simplified version - full implementation would use
    // sophisticated methods like:
    // - Graded meshes near singularities
    // - Coordinate transformations
    // - Subtraction of singularities

    // For now, use high-order Gauss quadrature with element subdivision
    GaussQuadrature gauss(8);

    // Subdivide element if it contains or is near the singularity
    Point2D midpoint = point1 * 0.5 + point2 * 0.5;

    Complex result1 = gauss.integrateLineElement(point1, midpoint, integrand);
    Complex result2 = gauss.integrateLineElement(midpoint, point2, integrand);

    return result1 + result2;
  }
};

/**
 * Integration strategy interface
 */
class IntegrationStrategy {
public:
  virtual ~IntegrationStrategy() = default;

  virtual Complex
  integrate(const Point2D &point1,
            const Point2D &point2,
            std::function<Complex(const Point2D &)> integrand) const = 0;

  virtual Complex
  integrateSingular(const Point2D &point1,
                    const Point2D &point2,
                    std::function<Complex(const Point2D &)> integrand,
                    const Point2D &singularity_point) const = 0;
};

/**
 * Standard Gauss quadrature strategy
 */
class StandardIntegration : public IntegrationStrategy {
public:
  explicit StandardIntegration(int order = 4) : gauss_(order) {
  }

  Complex
  integrate(const Point2D &point1,
            const Point2D &point2,
            std::function<Complex(const Point2D &)> integrand) const override {
    return gauss_.integrateLineElement(point1, point2, integrand);
  }

  Complex integrateSingular(const Point2D &point1,
                            const Point2D &point2,
                            std::function<Complex(const Point2D &)> integrand,
                            const Point2D &singularity_point) const override {
    AdaptiveIntegration adaptive;
    return adaptive.integrateAdaptive(point1, point2, integrand,
                                      singularity_point);
  }

private:
  GaussQuadrature gauss_;
};

} // namespace bem::domain
