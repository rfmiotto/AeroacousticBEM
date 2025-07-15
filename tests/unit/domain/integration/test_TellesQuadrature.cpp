#include "domain/integration/TellesQuadrature.hpp"
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

TEST(TellesQuadratureTest, CorrectNumberOfPointsPerOrder) {
  std::vector<std::pair<int, std::size_t>> order_points = {
      {2, 2}, {3, 3}, {4, 4}, {8, 8}};

  for (const auto &[order, expected_count] : order_points) {
    TellesQuadrature quad(order);
    EXPECT_EQ(quad.getPoints().size(), expected_count)
        << "Failed for order " << order;
  }
}

TEST(TellesQuadratureTest, ConstantFunction) {
  Element element = createUnitPanel();
  auto f = [](const Point2D &) -> Complex { return {1.0, 0.0}; };

  for (int order : {1, 2, 3, 4, 5, 6, 8}) {
    TellesQuadrature quad(order);
    Complex result = quad.integrate(element, element.midpoint(), f);
    EXPECT_NEAR(result.real(), 1.0, INTEGRATION_TOLERANCE)
        << "Failed on real part for order " << order;
    EXPECT_NEAR(result.imag(), 0.0, INTEGRATION_TOLERANCE)
        << "Failed on imag part for order " << order;
  }
}

TEST(TellesQuadratureTest, LinearFunction) {
  Element element = createUnitPanel();
  auto f = [](const Point2D &pt) -> Complex { return {pt.x, 0.0}; };

  for (int order : {1, 2, 3, 4, 5, 6, 8}) {
    TellesQuadrature quad(order);
    Complex result = quad.integrate(element, element.midpoint(), f);
    EXPECT_NEAR(result.real(), 0.5, INTEGRATION_TOLERANCE)
        << "Failed on real part for order " << order;
    EXPECT_NEAR(result.imag(), 0.0, INTEGRATION_TOLERANCE)
        << "Failed on imag part for order " << order;
  }
}

TEST(TellesQuadratureTest, ImaginaryFunction) {
  Element element = createUnitPanel();
  auto f = [](const Point2D &pt) -> Complex { return {0.0, pt.x}; };

  for (int order : {1, 2, 3, 4, 5, 6, 8}) {
    TellesQuadrature quad(order);
    Complex result = quad.integrate(element, element.midpoint(), f);
    EXPECT_NEAR(result.real(), 0.0, INTEGRATION_TOLERANCE)
        << "Failed on real part for order " << order;
    EXPECT_NEAR(result.imag(), 0.5, INTEGRATION_TOLERANCE)
        << "Failed on imag part for order " << order;
  }
}

TEST(TellesQuadratureTest, ThrowsForUnsupportedOrder) {
  EXPECT_THROW(TellesQuadrature(0), BEMIntegrationException);
  EXPECT_THROW(TellesQuadrature(7), BEMIntegrationException);
  EXPECT_THROW(TellesQuadrature(9), BEMIntegrationException);
}
