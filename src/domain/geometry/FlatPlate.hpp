#pragma once

#include "foundation/types/GeometryTypes.hpp"
#include <vector>

namespace bem::domain::geometry {

using bem::types::Panel;
using bem::types::PanelRegion;
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
  FlatPlate2D(std::size_t nAirfoilPanels,
              std::size_t nVerticalPanels,
              Real thickness,
              Real domainLength);

  [[nodiscard]] const std::vector<Point2D> &getPoints() const {
    return points_;
  }

  [[nodiscard]] const std::vector<Panel> &getPanels() const {
    return panels_;
  }

  [[nodiscard]] const std::vector<PanelRegion> &getTaggedPanels() const {
    return taggedPanels_;
  }

private:
  std::vector<Point2D> points_;
  std::vector<Panel> panels_;
  std::vector<PanelRegion> taggedPanels_;

  void generatePlate(std::size_t nAirfoilPanels,
                     std::size_t nVerticalPanels,
                     Real thickness,
                     Real domainLength);
};

} // namespace bem::domain::geometry
