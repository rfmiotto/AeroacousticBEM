#include "domain/integration/GaussQuadrature.hpp"
#include "foundation/utils/Constants.hpp"
#include <gtest/gtest.h>

using namespace bem::domain::integration;
using namespace bem::types;
using namespace bem::foundation::utils::Constants;

namespace {
Element createUnitPanel() {
  Point2D p1(0.0, 0.0);
  Point2D p2(1.0, 0.0);
  return {p1, p2};
}
} // namespace

TEST(GaussLegendreTableTest, CorrectNumberOfPointsPerOrder) {
  std::vector<int> orders = {1, 2, 3, 4, 5, 6, 8};
  for (int order : orders) {
    auto pts = GaussLegendreTable::getPoints(order);
    EXPECT_EQ(static_cast<std::size_t>(order), pts.size())
        << "Unexpected number of reference points for order " << order;
  }
}

// Test integral of f(x) = 1 over [0,1] → expect: 1.0
TEST(GaussQuadratureTest, ConstantFunction) {
  Element element = createUnitPanel();

  auto integrand = [](const QuadraturePoint & /*qp*/,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    val << Complex(1.0, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 4;
  GaussQuadrature quad(ip);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 1.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
}

// Test integral of f(x) = x over [0,1] → expect: 0.5
TEST(GaussQuadratureTest, LinearFunction) {
  Element element = createUnitPanel();

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    // qp.point is the physical coordinate of the quadrature point
    val << Complex(qp.point.x, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 6;
  GaussQuadrature quad(ip);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 0.5, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
}

// Test integral of f(x) = i*x over [0,1] → expect: 0.5i
TEST(GaussQuadratureTest, ImaginaryFunction) {
  Element element = createUnitPanel();

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    val << Complex(0.0, qp.point.x); // i * x
    return val;
  };

  IntegrationParameters ip;
  ip.order = 8;
  GaussQuadrature quad(ip);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 0.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.5, INTEGRATION_TOLERANCE);
}

// ∫_0^1 x^2 dx = 1/3
TEST(GaussQuadratureTest, QuadraticFunction) {
  Element element = createUnitPanel();

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    const double x = qp.point.x;
    val << Complex(x * x, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 6;
  GaussQuadrature quad(ip);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 1.0 / 3.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
}

// Vectorized test: return two entries simultaneously
TEST(GaussQuadratureTest, MultiDOFVectorizedIntegrand) {
  Element element = createUnitPanel();

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(2);
    const double x = qp.point.x;
    val << Complex(1.0, 0.0), Complex(x, 0.0); // two "shape functions"
    return val;
  };

  IntegrationParameters ip;
  ip.order = 6;
  GaussQuadrature quad(ip);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 2);
  // First entry: ∫_0^1 1 dx = 1
  EXPECT_NEAR(result(0).real(), 1.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
  // Second entry: ∫_0^1 x dx = 0.5
  EXPECT_NEAR(result(1).real(), 0.5, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(1).imag(), 0.0, INTEGRATION_TOLERANCE);
}
