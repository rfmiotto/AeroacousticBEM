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
  points_.reserve((nx + 1) + (nVerticalElements) + (nx) +
                  (nVerticalElements - 1)); // last point != first point

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

  // Define reference point indices for stretching (clockwise)
  topLeft_ = 0;
  topAirfoilLE_ = nUpDown;
  topAirfoilTE_ = topAirfoilLE_ + nAirfoilElements;
  topRight_ = topAirfoilTE_ + nUpDown;
  bottomRight_ = topRight_ + nVerticalElements;
  bottomAirfoilTE_ = bottomRight_ + nUpDown;
  bottomAirfoilLE_ = bottomAirfoilTE_ + nAirfoilElements;
  bottomLeft_ = bottomAirfoilLE_ + nUpDown;

  buildElements();
}

void FlatPlate2D::buildElements() {
  // Create segments and classify their region (UPSTREAM, AIRFOIL, DOWNSTREAM)
  elements_.clear();

  std::size_t nPoints = points_.size();

  elements_.reserve(nPoints);

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

    elements_.push_back(seg);
  }
}

void FlatPlate2D::setPoints(std::vector<Point2D> &&newPoints) {
  if (newPoints.size() != points_.size()) {
    throw std::invalid_argument("FlatPlate2D::setPoints: newPoints must have "
                                "the same size as the existing points.");
  }

  points_ = std::move(newPoints);
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getTopUpstreamPoints() const {
  return std::span<const Point2D>{&points_[topLeft_],
                                  topAirfoilLE_ - topLeft_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getTopAirfoilPoints() const {
  return std::span<const Point2D>{&points_[topAirfoilLE_],
                                  topAirfoilTE_ - topAirfoilLE_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getTopDownstreamPoints() const {
  return std::span<const Point2D>{&points_[topAirfoilTE_],
                                  topRight_ - topAirfoilTE_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getSideDownstreamPoints() const {
  return std::span<const Point2D>{&points_[topRight_],
                                  bottomRight_ - topRight_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getBottomDownstreamPoints() const {
  return std::span<const Point2D>{&points_[bottomRight_],
                                  bottomAirfoilTE_ - bottomRight_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getBottomAirfoilPoints() const {
  return std::span<const Point2D>{&points_[bottomAirfoilTE_],
                                  bottomAirfoilLE_ - bottomAirfoilTE_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getBottomUpstreamPoints() const {
  return std::span<const Point2D>{&points_[bottomAirfoilLE_],
                                  bottomLeft_ - bottomAirfoilLE_ + 1};
}

[[nodiscard]] std::span<const Point2D>
FlatPlate2D::getSideUpstreamPoints() const {
  return std::span<const Point2D>{&points_[bottomLeft_],
                                  points_.size() - bottomLeft_};
}
} // namespace bem::domain::geometry
