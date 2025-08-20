#include "domain/integration/TellesQuadrature.hpp"
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

// --- Core Tests ---------------------------------------------------------

TEST(TellesQuadratureTest, ConstantFunction) {
  Element element = createUnitPanel();
  Point2D x_star(0.5, 0.1); // near element

  auto integrand = [](const QuadraturePoint & /*qp*/,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    val << Complex(1.0, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 6;
  TellesQuadrature quad(ip, x_star);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 1.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
}

TEST(TellesQuadratureTest, LinearFunction) {
  Element element = createUnitPanel();
  Point2D x_star(0.25, 0.05); // shift collocation closer to edge

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    val << Complex(qp.point.x, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 8;
  TellesQuadrature quad(ip, x_star);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 0.5, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
}

TEST(TellesQuadratureTest, QuadraticFunction) {
  Element element = createUnitPanel();
  Point2D x_star(0.9, 0.2); // near endpoint

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    double x = qp.point.x;
    val << Complex(x * x, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 6;
  TellesQuadrature quad(ip, x_star);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 1.0 / 3.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.0, INTEGRATION_TOLERANCE);
}

TEST(TellesQuadratureTest, ImaginaryFunction) {
  Element element = createUnitPanel();
  Point2D x_star(0.5, 0.01); // almost on top of panel

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(1);
    val << Complex(0.0, qp.point.x); // i * x
    return val;
  };

  IntegrationParameters ip;
  ip.order = 6;
  TellesQuadrature quad(ip, x_star);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NEAR(result(0).real(), 0.0, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result(0).imag(), 0.5, INTEGRATION_TOLERANCE);
}

TEST(TellesQuadratureTest, MultiDOFIntegrand) {
  Element element = createUnitPanel();
  Point2D x_star(0.75, 0.15);

  auto integrand = [](const QuadraturePoint &qp,
                      const Element & /*e*/) -> Eigen::ArrayXcd {
    Eigen::ArrayXcd val(3);
    double x = qp.point.x;
    val << Complex(1.0, 0.0), Complex(x, 0.0), Complex(x * x, 0.0);
    return val;
  };

  IntegrationParameters ip;
  ip.order = 8;
  TellesQuadrature quad(ip, x_star);

  const Eigen::ArrayXcd result = quad.integrate(element, integrand);
  ASSERT_EQ(result.size(), 3);
  EXPECT_NEAR(result(0).real(), 1.0, INTEGRATION_TOLERANCE); // ∫1 dx = 1
  EXPECT_NEAR(result(1).real(), 0.5, INTEGRATION_TOLERANCE); // ∫x dx = 0.5
  EXPECT_NEAR(result(2).real(), 1.0 / 3.0,
              INTEGRATION_TOLERANCE); // ∫x² dx = 1/3
}
