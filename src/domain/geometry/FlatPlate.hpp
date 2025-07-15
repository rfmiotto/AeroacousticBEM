#pragma once

#include "foundation/types/GeometryTypes.hpp"
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

  [[nodiscard]] const std::vector<Point2D> &getPoints() const {
    return points_;
  }

  [[nodiscard]] const std::vector<Element> &getElements() const {
    return elements_;
  }

  [[nodiscard]] const std::vector<ElementRegion> &getTaggedElements() const {
    return taggedElements_;
  }

private:
  std::vector<Point2D> points_;
  std::vector<Element> elements_;
  std::vector<ElementRegion> taggedElements_;

  void generatePlate(std::size_t nAirfoilElements,
                     std::size_t nVerticalElements,
                     Real thickness,
                     Real domainLength);
};

} // namespace bem::domain::geometry
