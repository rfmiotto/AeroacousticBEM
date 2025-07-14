#include "domain/geometry/FlatPlate.hpp"
#include "foundation/utils/MathUtils.hpp"
#include "types/GeometryTypes.hpp"
#include <algorithm>
#include <gtest/gtest.h>
#include <ranges>

using namespace bem::types;
using namespace bem::domain::geometry;
using bem::foundation::utils::MathUtils::isEqual;

TEST(FlatPlate2DTest, GeneratesApproximateDomainDimensions) {
  std::size_t nAirfoil = 10;
  std::size_t nVertical = 4;
  Real thickness = 2.0;
  Real requestedLength = 5.0;

  FlatPlate2D plate(nAirfoil, nVertical, thickness, requestedLength);
  const auto &points = plate.getPoints();

  auto [xMinIt, xMaxIt] = std::ranges::minmax_element(points, {}, &Point2D::x);
  auto [yMinIt, yMaxIt] = std::ranges::minmax_element(points, {}, &Point2D::y);

  Real xMin = xMinIt->x;
  Real xMax = xMaxIt->x;
  Real yMin = yMinIt->y;
  Real yMax = yMaxIt->y;

  Real actualLength = xMax - xMin;
  Real dx =
      1.0 / static_cast<Real>(nAirfoil); // granularity of airfoil segments

  EXPECT_LT(xMin, 0.0) << "Domain should extend upstream of airfoil (x < 0)";
  EXPECT_GT(xMax, 1.0) << "Domain should extend downstream of airfoil (x > 1)";
  EXPECT_TRUE(isEqual(yMax - yMin, thickness)) << "Thickness mismatch";
  EXPECT_NEAR(actualLength, requestedLength, dx)
      << "Domain length is not approximately equal to the requested length";
}

TEST(FlatPlate2DTest, CorrectRegionClassification) {
  std::size_t nAirfoil = 20;
  std::size_t nVertical = 3;
  Real thickness = 1.5;
  Real domainLength = 7.0;

  FlatPlate2D plate(nAirfoil, nVertical, thickness, domainLength);
  const auto &segments = plate.getPanels();

  for (const auto &seg : segments) {
    Real mx = seg.midpoint().x;
    if (mx < 0.0) {
      EXPECT_EQ(seg.region, ElementRegion::UPSTREAM);
    } else if (mx <= 1.0) {
      EXPECT_EQ(seg.region, ElementRegion::AIRFOIL);
    } else {
      EXPECT_EQ(seg.region, ElementRegion::DOWNSTREAM);
    }
  }
}

TEST(FlatPlate2DTest, AirfoilSpansFromZeroToOne) {
  std::size_t nAirfoil = 16;
  std::size_t nVertical = 2;
  Real thickness = 1.0;
  Real domainLength = 4.0;

  FlatPlate2D plate(nAirfoil, nVertical, thickness, domainLength);
  const auto &segments = plate.getPanels();

  namespace views = std::views;

  auto is_airfoil = [](const auto &seg) {
    return seg.region == ElementRegion::AIRFOIL;
  };

  auto extract_x_coords = [](const auto &seg) {
    return std::array{seg.start.x, seg.end.x};
  };

  auto x_coords = segments | views::filter(is_airfoil) |
                  views::transform(extract_x_coords) | views::join;

  Real airfoilXMin = std::ranges::min(x_coords, std::less<>());
  Real airfoilXMax = std::ranges::max(x_coords, std::less<>());

  EXPECT_TRUE(isEqual(airfoilXMin, 0.0)) << "Airfoil must start at x = 0";
  EXPECT_TRUE(isEqual(airfoilXMax, 1.0)) << "Airfoil must end at x = 1";
}
