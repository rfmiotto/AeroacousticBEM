#pragma once

#include "foundation/types/GeometryTypes.hpp"
#include <span>
#include <vector>

namespace bem::domain::geometry {

using bem::types::Element;
using bem::types::ElementRegion;
using bem::types::Point2D;
using bem::types::Real;

/**
 * @brief Represents a 2D flat plate in a rectangular domain for boundary
 * element methods.
 *
 * The domain is divided into three logical regions:
 * - Upstream:   x ∈ [xBeg, 0)
 * - Airfoil:    x ∈ [0, 1]
 * - Downstream: x ∈ (1, xEnd]
 *
 * The plate is flat and lies along y = 0, with optional vertical
 * discretization.
 */
class FlatPlate2D {
public:
  FlatPlate2D(std::size_t nAirfoilElements,
              std::size_t nVerticalElements,
              Real thickness,
              Real domainLength);

  [[nodiscard]] std::vector<Point2D> &accessPoints() {
    return points_;
  }

  [[nodiscard]] const std::vector<Point2D> &getPoints() const {
    return points_;
  }

  [[nodiscard]] const std::vector<Element> &getElements() const {
    return elements_;
  }

  void setPoints(std::vector<Point2D> &&newPoints);

  void buildElements();

  [[nodiscard]] std::span<const Point2D> getTopUpstreamPoints() const;
  [[nodiscard]] std::span<const Point2D> getTopAirfoilPoints() const;
  [[nodiscard]] std::span<const Point2D> getTopDownstreamPoints() const;

  [[nodiscard]] std::span<const Point2D> getSideUpstreamPoints() const;
  [[nodiscard]] std::span<const Point2D> getSideDownstreamPoints() const;

  [[nodiscard]] std::span<const Point2D> getBottomUpstreamPoints() const;
  [[nodiscard]] std::span<const Point2D> getBottomAirfoilPoints() const;
  [[nodiscard]] std::span<const Point2D> getBottomDownstreamPoints() const;

private:
  std::vector<Point2D> points_;
  std::vector<Element> elements_;

  // Reference point indices
  std::size_t topLeft_{0};
  std::size_t topAirfoilLE_{0};
  std::size_t topAirfoilTE_{0};
  std::size_t topRight_{0};
  std::size_t bottomRight_{0};
  std::size_t bottomAirfoilTE_{0};
  std::size_t bottomAirfoilLE_{0};
  std::size_t bottomLeft_{0};

  void generatePlate(std::size_t nAirfoilElements,
                     std::size_t nVerticalElements,
                     Real thickness,
                     Real domainLength);
};

} // namespace bem::domain::geometry
