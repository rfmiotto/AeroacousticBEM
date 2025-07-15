#include "domain/integration/GaussQuadrature.hpp"
#include "foundation/exceptions/BEMException.hpp"
#include "foundation/utils/Constants.hpp"
#include <gtest/gtest.h>

using namespace bem::domain::integration;
using namespace bem::types;
using namespace bem::foundation::utils::Constants;
using namespace bem::foundation::exceptions;

namespace {
Element createUnitPanel() {
  Point2D p1(0.0, 0.0);
  Point2D p2(1.0, 0.0);
  return {p1, p2};
}
} // namespace

TEST(GaussQuadratureTest, CorrectNumberOfPointsPerOrder) {
  std::vector<std::pair<int, std::size_t>> order_points = {
      {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {8, 8}};

  for (const auto &[order, expected_count] : order_points) {
    GaussQuadrature quad(order);
    EXPECT_EQ(quad.getPoints().size(), expected_count)
        << "Failed for order " << order;
  }
}

// Test integral of f(x) = 1 over [0,1] → expect: 1.0
TEST(GaussQuadratureTest, ConstantFunction) {
  Element element = createUnitPanel();

  auto f = [](const Point2D &) -> Complex { return {1.0, 0.0}; };

  for (int order : {1, 2, 3, 4, 5, 6, 8}) {
    GaussQuadrature quad(order);
    Complex result = quad.integrate(element, element.midpoint(), f);
    EXPECT_NEAR(result.real(), 1.0, INTEGRATION_TOLERANCE)
        << "Failed on real part for order " << order;
    EXPECT_NEAR(result.imag(), 0.0, INTEGRATION_TOLERANCE)
        << "Failed on imag part for order " << order;
  }
}

// Test integral of f(x) = x over [0,1] → expect: 0.5
TEST(GaussQuadratureTest, LinearFunction) {
  Element element = createUnitPanel();

  auto f = [](const Point2D &pt) -> Complex { return {pt.x, 0.0}; };

  for (int order : {1, 2, 3, 4, 5, 6, 8}) {
    GaussQuadrature quad(order);
    Complex result = quad.integrate(element, element.midpoint(), f);
    EXPECT_NEAR(result.real(), 0.5, INTEGRATION_TOLERANCE)
        << "Failed on real part for order " << order;
    EXPECT_NEAR(result.imag(), 0.0, INTEGRATION_TOLERANCE)
        << "Failed on imag part for order " << order;
  }
}

// Test integral of f(x) = i*x over [0,1] → expect: 0.5i
TEST(GaussQuadratureTest, ImaginaryFunction) {
  Element element = createUnitPanel();

  auto f = [](const Point2D &pt) -> Complex { return {0.0, pt.x}; };

  for (int order : {1, 2, 3, 4, 5, 6, 8}) {
    GaussQuadrature quad(order);
    Complex result = quad.integrate(element, element.midpoint(), f);
    EXPECT_NEAR(result.real(), 0.0, INTEGRATION_TOLERANCE)
        << "Failed on real part for order " << order;
    EXPECT_NEAR(result.imag(), 0.5, INTEGRATION_TOLERANCE)
        << "Failed on imag part for order " << order;
  }
}

TEST(GaussQuadratureTest, ThrowsForUnsupportedOrder) {
  EXPECT_THROW(GaussQuadrature(7), BEMIntegrationException);
  EXPECT_THROW(GaussQuadrature(0), BEMIntegrationException);
  EXPECT_THROW(GaussQuadrature(9), BEMIntegrationException);
}
