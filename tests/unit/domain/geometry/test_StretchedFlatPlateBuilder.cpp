#include "domain/geometry/StretchedFlatPlateBuilder.hpp"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

using bem::domain::geometry::FlatPlate2D;
using bem::domain::geometry::makeStretchedPlate;
using bem::domain::geometry::StretchingParams;
using bem::types::Point2D;
using bem::types::Real;

namespace {

void writePointsToCSV(const std::vector<Point2D> &points,
                      const std::string &filename) {
  std::ofstream file(filename);
  if (!file) {
    return;
  }
  file << "x,y\n";
  for (const auto &pt : points) {
    file << pt.x << "," << pt.y << "\n";
  }
  file.close();
}

bool isIncreasing(const std::vector<double> &deltas) {
  for (std::size_t i = 1; i < deltas.size(); ++i) {
    if (std::abs(deltas[i]) < std::abs(deltas[i - 1]) - 1e-12) {
      return false;
    }
  }
  return true;
}

bool isDecreasing(const std::vector<double> &deltas) {
  for (std::size_t i = 1; i < deltas.size(); ++i) {
    if (std::abs(deltas[i]) > std::abs(deltas[i - 1]) + 1e-12) {
      return false;
    }
  }
  return true;
}

std::vector<double> computeDeltas(const std::vector<Point2D> &points) {
  std::vector<double> deltas;
  for (std::size_t i = 1; i < points.size(); ++i) {
    deltas.push_back(points[i].x - points[i - 1].x);
  }
  return deltas;
}

} // namespace

TEST(StretchedFlatPlateBuilderTest, AppliesCorrectStretchingNearEdges) {
  const std::size_t nAirfoil = 10;
  const std::size_t nVertical = 2;
  const Real thickness = 0.05;
  const Real domainLength = 3.0;

  StretchingParams params{.paramP = 1.9, .paramQ = 2.0};

  FlatPlate2D plate =
      makeStretchedPlate(nAirfoil, nVertical, thickness, domainLength, params);
  const auto &points = plate.getPoints();

  // Export all points to CSV for visual inspection
  writePointsToCSV(points, "stretched_plate_points.csv");

  // Get all regions
  const auto topUpstream = plate.getTopUpstreamPoints();
  const auto botUpstream = plate.getBottomUpstreamPoints();
  const auto topDownstream = plate.getTopDownstreamPoints();
  const auto botDownstream = plate.getBottomDownstreamPoints();

  const auto topAirfoil = plate.getTopAirfoilPoints();
  const auto botAirfoil = plate.getBottomAirfoilPoints();

  // Upstream
  auto deltaTopUp = computeDeltas({topUpstream.begin(), topUpstream.end()});
  auto deltaBotUp = computeDeltas({botUpstream.begin(), botUpstream.end()});

  EXPECT_TRUE(isDecreasing(deltaTopUp)) << "Top upstream not decreasing";
  EXPECT_TRUE(isIncreasing(deltaBotUp)) << "Bottom upstream not increasing";

  // Downstream
  auto deltaTopDown =
      computeDeltas({topDownstream.begin(), topDownstream.end()});
  auto deltaBotDown =
      computeDeltas({botDownstream.begin(), botDownstream.end()});

  EXPECT_TRUE(isIncreasing(deltaTopDown)) << "Top downstream not increasing";
  EXPECT_TRUE(isDecreasing(deltaBotDown)) << "Bottom downstream not decreasing";

  // Airfoil - symmetrically split
  auto splitSpanInHalfSymmetric = [](std::span<const Point2D> span) {
    const auto half = static_cast<std::ptrdiff_t>(span.size() / 2);
    return std::pair{std::vector<Point2D>{span.begin(), span.begin() + half},
                     std::vector<Point2D>{span.end() - half, span.end()}};
  };

  auto [topLeading, topTrailing] = splitSpanInHalfSymmetric(topAirfoil);
  auto [botTrailing, botLeading] = splitSpanInHalfSymmetric(botAirfoil);

  auto deltaTopLead = computeDeltas(topLeading);
  auto deltaTopTrail = computeDeltas(topTrailing);
  auto deltaBotLead = computeDeltas(botLeading);
  auto deltaBotTrail = computeDeltas(botTrailing);

  EXPECT_TRUE(isIncreasing(deltaTopLead)) << "Top leading edge not increasing";
  EXPECT_TRUE(isDecreasing(deltaTopTrail))
      << "Top trailing edge not decreasing";
  EXPECT_TRUE(isDecreasing(deltaBotLead))
      << "Bottom leading edge not decreasing";
  EXPECT_TRUE(isIncreasing(deltaBotTrail))
      << "Bottom trailing edge not increasing";

  // Sanity checks
  EXPECT_GT(points.size(), 0);

  // Verify that plate starts and ends at correct bounds
  const Real dx = 1.0 / static_cast<Real>(nAirfoil);
  const auto nUpDown =
      static_cast<std::size_t>(std::round((domainLength - 1.0) / (2.0 * dx)));
  Real expectedBeg = -dx * static_cast<Real>(nUpDown);
  Real expectedEnd = 1.0 + (dx * static_cast<Real>(nUpDown));

  EXPECT_NEAR(points.front().x, expectedBeg, 1e-8);
  EXPECT_NEAR(topDownstream.back().x, expectedEnd, 1e-8);
}
