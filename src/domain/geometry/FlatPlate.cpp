#include "domain/geometry/FlatPlate.hpp"
#include <cstdlib>

namespace bem::domain::geometry {
using namespace bem::types;

FlatPlate2D::FlatPlate2D(std::size_t nAirfoilElements,
                         std::size_t nVerticalElements,
                         Real thickness,
                         Real domainLength) {
  generatePlate(nAirfoilElements, nVerticalElements, thickness, domainLength);
}
void FlatPlate2D::generatePlate(std::size_t nAirfoilElements,
                                std::size_t nVerticalElements,
                                Real thickness,
                                Real domainLength) {
  if (domainLength <= 1.0) {
    throw std::invalid_argument("Domain length must be > 1.0 to include "
                                "airfoil + upstream/downstream.");
  }
  if (nAirfoilElements < 1) {
    throw std::invalid_argument("nAirfoil must be >= 1.");
  }
  if (thickness <= 0.0) {
    throw std::invalid_argument("Thickness must be positive.");
  }
  if (nVerticalElements < 1) {
    throw std::invalid_argument("nVerticalSegments must be >= 1.");
  }

  Real dx = 1.0 / static_cast<Real>(nAirfoilElements);
  auto nUpDown =
      static_cast<std::size_t>(std::round((domainLength - 1.0) / (2.0 * dx)));

  std::size_t nx = (2 * nUpDown) + nAirfoilElements;
  Real xBeg = -dx * static_cast<Real>(nUpDown);
  Real xEnd = 1.0 + (dx * static_cast<Real>(nUpDown));
  Real yBeg = -0.5 * thickness;
  Real yEnd = 0.5 * thickness;
  Real dy = (yEnd - yBeg) / static_cast<Real>(nVerticalElements);

  points_.clear();
  panels_.clear();
  taggedPanels_.clear();

  // 1. Top side: left → right
  for (std::size_t i = 0; i <= nx; ++i) {
    Real x = xBeg + (static_cast<Real>(i) * dx);
    points_.emplace_back(x, yEnd);
  }

  // 2. Right side: top → bottom
  for (std::size_t j = 1; j <= nVerticalElements; ++j) {
    Real y = yEnd - (static_cast<Real>(j) * dy);
    points_.emplace_back(xEnd, y);
  }

  // 3. Bottom side: right → left
  for (std::size_t i = 1; i <= nx; ++i) {
    Real x = xEnd - (static_cast<Real>(i) * dx);
    points_.emplace_back(x, yBeg);
  }

  // 4. Left side: bottom → top (excluding top-left corner)
  for (std::size_t j = 1; j < nVerticalElements; ++j) {
    Real y = yBeg + (static_cast<Real>(j) * dy);
    points_.emplace_back(xBeg, y);
  }

  // Create segments and classify their region (UPSTREAM, AIRFOIL, DOWNSTREAM)
  std::size_t nPoints = points_.size();
  for (std::size_t i = 0; i < nPoints; ++i) {
    const auto &p1 = points_[i];
    const auto &p2 = points_[(i + 1) % nPoints]; // closed loop
    Element seg(p1, p2);

    // Determine region
    Real mx = seg.midpoint().x;
    if (mx < 0.0) {
      seg.region = ElementRegion::UPSTREAM;
    } else if (mx <= 1.0) {
      seg.region = ElementRegion::AIRFOIL;
    } else {
      seg.region = ElementRegion::DOWNSTREAM;
    }

    panels_.push_back(seg);
  }
}
} // namespace bem::domain::geometry
