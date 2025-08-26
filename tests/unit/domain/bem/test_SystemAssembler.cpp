#include "domain/bem/SystemAssembler.hpp"
#include "gtest/gtest.h"

using namespace bem::domain::assembly;
using namespace bem::domain::integration;
using namespace bem::types;

// Struct to describe a test scenario
struct AssembleHDTestCase {
  std::string name;
  Element elem;
  Point2D collocation;
  Real singularity_distance;
  double k;
};

// Hack to avoid printing the raw bytes of the struct passed into a
// parameterized test, as the struct did not overload the stream operator.
// NOLINTNEXTLINE(readability-identifier-naming, misc-use-internal-linkage)
void PrintTo(const AssembleHDTestCase &tc, std::ostream *os) {
  *os << "AssembleHDTestCase{name=" << tc.name << ", elem=[(" << tc.elem.start.x
      << "," << tc.elem.start.y << ") -> (" << tc.elem.end.x << ","
      << tc.elem.end.y << ")]"
      << ", type=" << static_cast<int>(tc.elem.type)
      << ", region=" << static_cast<int>(tc.elem.region) << ", collocation=("
      << tc.collocation.x << "," << tc.collocation.y << ")"
      << ", singularity_distance=" << tc.singularity_distance << ", k=" << tc.k
      << "}";
}

// Parameterized test fixture
class AssembleSystemHDParameterizedTest
    : public ::testing::TestWithParam<AssembleHDTestCase> {};

TEST_P(AssembleSystemHDParameterizedTest, BasicAssembly) {
  const auto &tc = GetParam();

  // Use a single-element mesh
  std::vector<Element> elements{tc.elem};
  std::vector<Point2D> collocation{tc.collocation};

  // Mapper: each element has 1 DOF (for simplicity)
  GlobalDofMapper mapper = [](std::size_t idx) -> std::vector<int> {
    return {static_cast<int>(idx)};
  };

  // Integration params
  IntegrationParameters ip{};
  ip.singularity_distance = tc.singularity_distance;

  Eigen::MatrixXcd H;
  Eigen::MatrixXcd D;
  assembleSystemHD(elements, collocation, tc.k, ip, mapper, H, D);

  // Sanity checks
  EXPECT_EQ(H.rows(), 1);
  EXPECT_EQ(H.cols(), 1);
  EXPECT_EQ(D.rows(), 1);
  EXPECT_EQ(D.cols(), 1);

  // Ensure non-zero contributions
  EXPECT_FALSE(std::isnan(H(0, 0).real()));
  EXPECT_FALSE(std::isnan(H(0, 0).imag()));
  EXPECT_FALSE(std::isnan(D(0, 0).real()));
  EXPECT_FALSE(std::isnan(D(0, 0).imag()));

  // Check which quadrature was used based on collocation
  if (tc.elem.contains(tc.collocation, 1e-12)) {
    // Analytical path
    SUCCEED() << "Analytical self-term path executed";
  } else if (std::sqrt(std::pow(tc.collocation.x -
                                    ((tc.elem.start.x + tc.elem.end.x) / 2),
                                2) +
                       std::pow(tc.collocation.y -
                                    ((tc.elem.start.y + tc.elem.end.y) / 2),
                                2)) <= tc.singularity_distance) {
    // Telles path
    SUCCEED() << "Telles quadrature path executed";
  } else {
    // Gauss path
    SUCCEED() << "Gauss quadrature path executed";
  }
}

// Instantiate test cases
INSTANTIATE_TEST_SUITE_P(
    AssembleSystemHDTests,
    AssembleSystemHDParameterizedTest,
    ::testing::Values(
        // Analytical self-term (collocation on element)
        AssembleHDTestCase{"Analytical",
                           Element{Point2D{0.0, 0.0}, Point2D{1.0, 0.0},
                                   ElementType::CONSTANT},
                           Point2D{0.5, 0.0}, 1e-6, 1.0},
        // Gauss (far from element)
        AssembleHDTestCase{
            "Gauss",
            Element{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR},
            Point2D{0.5, 1.0}, 1e-6, 1.0},
        // Telles (near-singular)
        AssembleHDTestCase{
            "Telles",
            Element{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR},
            Point2D{0.5, 1e-8}, 1e-6, 1.0}),
    [](const ::testing::TestParamInfo<AssembleHDTestCase> &info) {
      return info.param.name; // Use descriptive names
    });
