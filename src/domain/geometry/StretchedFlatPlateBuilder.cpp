#include "StretchedFlatPlateBuilder.hpp"

namespace bem::domain::geometry {

FlatPlate2D makeStretchedPlate(std::size_t nAirfoilElements,
                               std::size_t nVerticalElements,
                               Real thickness,
                               Real domainLength,
                               const StretchingParams &params) {
  // 1. Create flat plate with uniform point distribution
  FlatPlate2D plate(nAirfoilElements, nVerticalElements, thickness,
                    domainLength);
  std::vector<Point2D> &stretchedPoints = plate.accessPoints();

  // Split span into equal halves (symmetrical), discard middle if odd
  auto splitSpanInHalfSymmetric = [](std::span<const Point2D> span) {
    const std::size_t half = span.size() / 2;
    return std::pair{span.subspan(0, half),
                     span.subspan(span.size() - half, half)};
  };

  // Helper to apply stretching to a span of points in-place
  auto stretchSpan = [&](std::span<const Point2D> span, Real xBeg, Real xEnd,
                         bool concentrateAtEnd) {
    const int n = static_cast<int>(span.size());
    std::ptrdiff_t offset = span.data() - stretchedPoints.data();

    StretchingParams localParams = params;
    localParams.concentrateAtEnd = concentrateAtEnd;

    StretchingFunction stretch(xBeg, xEnd, n, localParams);
    const auto nodes = stretch.generateNodes();

    for (int i = 0; i < n; ++i) {
      stretchedPoints[offset + i].x = nodes[i];
    }
  };

  // 2. Apply stretching to each region

  // -- UPSTREAM
  {
    auto topUpstream = plate.getTopUpstreamPoints();
    auto botUpstream = plate.getBottomUpstreamPoints();

    stretchSpan(topUpstream, topUpstream.front().x, topUpstream.back().x, true);
    stretchSpan(botUpstream, botUpstream.front().x, botUpstream.back().x,
                false);
  }

  // -- AIRFOIL
  {
    auto [topLeft, topRight] =
        splitSpanInHalfSymmetric(plate.getTopAirfoilPoints());
    auto [botRight, botLeft] =
        splitSpanInHalfSymmetric(plate.getBottomAirfoilPoints());

    stretchSpan(topLeft, topLeft.front().x, topLeft.back().x, false);
    stretchSpan(topRight, topRight.front().x, topRight.back().x, true);
    stretchSpan(botRight, botRight.front().x, botRight.back().x, false);
    stretchSpan(botLeft, botLeft.front().x, botLeft.back().x, true);
  }

  // -- DOWNSTREAM
  {
    auto topDownstream = plate.getTopDownstreamPoints();
    auto botDownstream = plate.getBottomDownstreamPoints();

    stretchSpan(topDownstream, topDownstream.front().x, topDownstream.back().x,
                false);
    stretchSpan(botDownstream, botDownstream.front().x, botDownstream.back().x,
                true);
  }

  return plate;
}

} // namespace bem::domain::geometry
