#include "domain/integration/AnalyticalConstantSelfTerm.hpp"
#include "foundation/utils/Constants.hpp"
#include <complex>
#include <gtest/gtest.h>

using namespace bem::domain::integration;
using namespace bem::foundation::utils;
using namespace bem::types;

namespace {

TEST(AnalyticalConstantQuadratureTest, SingleLayerAndDoubleLayerConsistency) {
  // Define a simple element: segment from (0,0) to (2,0)
  Element elem{Point2D{0.0, 0.0}, Point2D{2.0, 0.0}};

  // Length should be 2.0
  double length = elem.length();
  ASSERT_NEAR(length, 2.0, 1e-14);

  // Define wavenumber
  double k = 1.0;
  AnalyticalConstantQuadrature quad(k);

  // Call integrate() with dummy integrand (not used by
  // AnalyticalConstantQuadrature)
  auto result =
      quad.integrate(elem, [](const QuadraturePoint &, const Element &) {
        return Eigen::ArrayXcd::Zero(1);
      });

  // Extract results
  const auto &single = quad.singleLayer();
  const auto &dbl = quad.doubleLayer();

  ASSERT_EQ(single.size(), 1);
  ASSERT_EQ(dbl.size(), 1);

  // Expected values
  double term =
      (1.0 - Constants::EULER_MASCHERONI - std::log(k * length * 0.25)) *
      Constants::INV_2PI;
  std::complex<double> expected_single(term * length, 0.25 * length);
  std::complex<double> expected_double(0.5, 0.0);

  // Compare with tolerances
  EXPECT_NEAR(single(0).real(), expected_single.real(), 1e-12);
  EXPECT_NEAR(single(0).imag(), expected_single.imag(), 1e-12);
  EXPECT_NEAR(dbl(0).real(), expected_double.real(), 1e-12);
  EXPECT_NEAR(dbl(0).imag(), expected_double.imag(), 1e-12);
}

TEST(AnalyticalConstantQuadratureTest, ZeroLengthElement) {
  // Degenerate element: zero length
  Element elem{Point2D{0.0, 0.0}, Point2D{0.0, 0.0}};

  double k = 2.0;
  AnalyticalConstantQuadrature quad(k);

  auto result =
      quad.integrate(elem, [](const QuadraturePoint &, const Element &) {
        return Eigen::ArrayXcd::Zero(1);
      });

  const auto &single = quad.singleLayer();
  const auto &dbl = quad.doubleLayer();

  // Zero length => both single- and double-layer should collapse to (0,0)
  EXPECT_NEAR(single(0).real(), 0.0, 1e-14);
  EXPECT_NEAR(single(0).imag(), 0.0, 1e-14);

  EXPECT_NEAR(dbl(0).real(), 0.0, 1e-14);
  EXPECT_NEAR(dbl(0).imag(), 0.0, 1e-14);
}

} // namespace
