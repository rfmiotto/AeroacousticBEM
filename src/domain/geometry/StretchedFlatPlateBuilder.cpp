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
  std::vector<Point2D> stretchedPoints = plate.getPoints();

  // Helper to apply stretching in-place to a span of points
  auto applyStretching = [](std::span<const Point2D> pointsToStretch,
                            std::vector<Point2D> &allPoints, Real xBeg,
                            Real xEnd, bool reverse,
                            const StretchingParams &params) {
    const int n = static_cast<int>(pointsToStretch.size());
    const std::ptrdiff_t offset = pointsToStretch.data() - allPoints.data();

    StretchingFunction stretch(xBeg, xEnd, n, params);
    const auto nodes = stretch.generateNodes();

    for (int i = 0; i < n; ++i) {
      const std::size_t idx = offset + i;
      allPoints[idx].x = nodes[reverse ? (n - 1 - i) : i];
    }
  };

  // Helper to split a span into two halves
  auto splitSpanInHalf = [](std::span<const Point2D> span) {
    const std::size_t half = span.size() / 2;
    return std::pair{span.first(half), span.last(span.size() - half)};
  };

  // 2. Apply stretching to each region

  // -- UPSTREAM
  {
    auto topUpstream = plate.getTopUpstreamPoints();
    auto idx = topUpstream.data() - stretchedPoints.data();
    applyStretching(topUpstream, stretchedPoints, stretchedPoints[idx].x, 0.0,
                    false, params);

    auto bottomUpstream = plate.getBottomUpstreamPoints();
    idx = bottomUpstream.data() - stretchedPoints.data();
    applyStretching(bottomUpstream, stretchedPoints, stretchedPoints[idx].x,
                    0.0, true, params);
  }

  // -- AIRFOIL
  {
    auto [topLeft, topRight] = splitSpanInHalf(plate.getTopAirfoilPoints());
    auto [botRight, botLeft] =
        splitSpanInHalf(plate.getBottomAirfoilPoints()); // mirrored

    applyStretching(topLeft, stretchedPoints, 0.0, 0.5, false, params);
    applyStretching(topRight, stretchedPoints, 0.5, 1.0, false, params);
    applyStretching(botLeft, stretchedPoints, 0.0, 0.5, true, params);
    applyStretching(botRight, stretchedPoints, 0.5, 1.0, true, params);
  }

  // -- DOWNSTREAM
  {
    auto topDownstream = plate.getTopDownstreamPoints();
    auto idx = topDownstream.data() - stretchedPoints.data();
    applyStretching(topDownstream, stretchedPoints, 1.0, stretchedPoints[idx].x,
                    false, params);

    auto bottomDownstream = plate.getBottomDownstreamPoints();
    idx = bottomDownstream.data() - stretchedPoints.data();
    applyStretching(bottomDownstream, stretchedPoints, 1.0,
                    stretchedPoints[idx].x, true, params);
  }

  // 3. Update plate with modified points
  plate.setPoints(std::move(stretchedPoints));
  return plate;
}

} // namespace bem::domain::geometry
