#include "domain/integration/ElementKernelContrib.hpp"
#include "domain/physics/GreensFunctions.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <ostream>

using namespace bem::domain::integration;
using namespace bem::types;

// ------------------ Basic Functionality Test ------------------
TEST(ComputeVectorizedContribAtQPTest, BasicFunctionality) {
  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}};
  QuadraturePoint qp{0.0, 1.0, Point2D{0.5, 0.0}};
  Point2D collocation{0.5, 0.5};
  double k = 2.0;

  const Eigen::Index nDOF =
      static_cast<Eigen::Index>(elem.shapeFunction(qp.xi).size());

  // SINGLE_LAYER_G
  auto vecG = computeVectorizedContribAtQP(qp, elem, collocation, k,
                                           KernelType::SINGLE_LAYER_G);
  EXPECT_EQ(vecG.size(), nDOF);
  for (Eigen::Index i = 0; i < nDOF; ++i) {
    EXPECT_FALSE(std::isnan(vecG(i).real()));
    EXPECT_FALSE(std::isnan(vecG(i).imag()));
  }

  // DOUBLE_LAYER_DGDN
  auto vecDG = computeVectorizedContribAtQP(qp, elem, collocation, k,
                                            KernelType::DOUBLE_LAYER_DGDN);
  EXPECT_EQ(vecDG.size(), nDOF);
  for (Eigen::Index i = 0; i < nDOF; ++i) {
    EXPECT_FALSE(std::isnan(vecDG(i).real()));
    EXPECT_FALSE(std::isnan(vecDG(i).imag()));
  }

  // Ensure the single- and double-layer contributions differ
  bool all_equal = true;
  for (Eigen::Index i = 0; i < nDOF; ++i) {
    if (vecG(i) != vecDG(i)) {
      all_equal = false;
      break;
    }
  }
  EXPECT_FALSE(all_equal);
}

// ------------------ Analytical Comparison ------------------
struct GreenTestParams {
  Element elem;
  QuadraturePoint qp;
  Point2D collocation;
  double k;
};

// Hack to avoid printing the raw bytes of the struct passed into a
// parameterized test, as the struct did not overload the stream operator.
// NOLINTNEXTLINE(readability-identifier-naming, misc-use-internal-linkage)
inline void PrintTo(const GreenTestParams &p, ::std::ostream *os) {
  *os << "Element[(" << p.elem.start.x << "," << p.elem.start.y << ") -> ("
      << p.elem.end.x << "," << p.elem.end.y << ")] "
      << "QP(" << p.qp.point.x << "," << p.qp.point.y << ") "
      << "Collocation(" << p.collocation.x << "," << p.collocation.y << ") "
      << "k=" << p.k;
}

class ComputeVectorizedContribAtQPAnalyticalTest
    : public ::testing::TestWithParam<GreenTestParams> {};

TEST_P(ComputeVectorizedContribAtQPAnalyticalTest,
       SingleLayerMatchesAnalytical) {
  const auto &p = GetParam();
  auto result = computeVectorizedContribAtQP(p.qp, p.elem, p.collocation, p.k,
                                             KernelType::SINGLE_LAYER_G);

  auto N = p.elem.shapeFunction(p.qp.xi);
  for (Eigen::Index i = 0; i < result.size(); ++i) {
    std::complex<double> expected =
        bem::domain::physics::greensFunction2D(p.qp.point, p.collocation, p.k) *
        N[i];
    EXPECT_NEAR(result(i).real(), expected.real(), 1e-12);
    EXPECT_NEAR(result(i).imag(), expected.imag(), 1e-12);
  }
}

TEST_P(ComputeVectorizedContribAtQPAnalyticalTest,
       DoubleLayerMatchesAnalytical) {
  const auto &p = GetParam();
  auto result = computeVectorizedContribAtQP(p.qp, p.elem, p.collocation, p.k,
                                             KernelType::DOUBLE_LAYER_DGDN);

  auto N = p.elem.shapeFunction(p.qp.xi);
  Vector2D n = p.elem.normal();
  for (Eigen::Index i = 0; i < result.size(); ++i) {
    std::complex<double> expected =
        bem::domain::physics::greensFunctionNormalDerivative2D(
            p.qp.point, p.collocation, n, p.k) *
        N[i];
    EXPECT_NEAR(result(i).real(), expected.real(), 1e-12);
    EXPECT_NEAR(result(i).imag(), expected.imag(), 1e-12);
  }
}

// Parameterized test cases
INSTANTIATE_TEST_SUITE_P(
    GreenFunctionTests,
    ComputeVectorizedContribAtQPAnalyticalTest,
    ::testing::Values(
        // Non-singular
        GreenTestParams{Element{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}},
                        QuadraturePoint{0.0, 1.0, Point2D{0.5, 0.0}},
                        Point2D{0.5, 0.5}, 2.0},
        // Near-singular (collocation very close to element)
        GreenTestParams{Element{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}},
                        QuadraturePoint{0.0, 1.0, Point2D{0.5, 0.0}},
                        Point2D{0.5, 1e-8}, 2.0}));
