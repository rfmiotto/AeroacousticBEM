#include "domain/geometry/Stretching.hpp"
#include "utils/Constants.hpp"
#include <gtest/gtest.h>

using namespace bem::domain::geometry;
using namespace bem::foundation::utils::Constants;

TEST(StretchingFunctionTest, GeneratesCorrectNumberOfNodes) {
  StretchingParams params{.paramP = 0.5, .paramQ = 3.0};
  StretchingFunction stretch(0.0, 1.0, 11, params);
  std::vector<double> nodes = stretch.generateNodes();

  ASSERT_EQ(nodes.size(), 11);
}

TEST(StretchingFunctionTest,
     NodesRespectIncreasingDomainBoundsAndConcentrationAtEnd) {
  // --- Default stretching
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = true};
  double xBeg = 0.0;
  double xEnd = 1.0;
  int numPoints = 11;

  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto nodes = stretch.generateNodes();

  ASSERT_NEAR(nodes.front(), xBeg, GEOMETRY_TOLERANCE);
  ASSERT_NEAR(nodes.back(), xEnd, GEOMETRY_TOLERANCE);
}

TEST(StretchingFunctionTest,
     NodesRespectIncreasingDomainBoundsAndConcentrationAtBeg) {
  // --- Default stretching
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = false};
  double xBeg = 0.0;
  double xEnd = 1.0;
  int numPoints = 11;

  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto nodes = stretch.generateNodes();

  ASSERT_NEAR(nodes.front(), xBeg, GEOMETRY_TOLERANCE);
  ASSERT_NEAR(nodes.back(), xEnd, GEOMETRY_TOLERANCE);
}

TEST(StretchingFunctionTest,
     NodesRespectDecreasingDomainBoundsAndConcentrationAtEnd) {
  // --- Default stretching
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = true};
  double xBeg = 1.0;
  double xEnd = 0.0;
  int numPoints = 11;

  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto nodes = stretch.generateNodes();

  ASSERT_NEAR(nodes.front(), xBeg, GEOMETRY_TOLERANCE);
  ASSERT_NEAR(nodes.back(), xEnd, GEOMETRY_TOLERANCE);
}

TEST(StretchingFunctionTest,
     NodesRespectDecreasingDomainBoundsAndConcentrationAtBeg) {
  // --- Default stretching
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = false};
  double xBeg = 1.0;
  double xEnd = 0.0;
  int numPoints = 11;

  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto nodes = stretch.generateNodes();

  ASSERT_NEAR(nodes.front(), xBeg, GEOMETRY_TOLERANCE);
  ASSERT_NEAR(nodes.back(), xEnd, GEOMETRY_TOLERANCE);
}

TEST(StretchingFunctionTest, ConcentratesAtBeginningWhenSpecified) {
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = false};
  StretchingFunction stretch(0.0, 1.0, 20, params);
  auto nodes = stretch.generateNodes();

  // Expect nodes to be closer to each other at the beginning (0.0)
  double first_gap = nodes[1] - nodes[0];
  double last_gap = nodes[nodes.size() - 1] - nodes[nodes.size() - 2];
  ASSERT_LT(first_gap, last_gap)
      << "Spacing at beginning is not smaller than at end as expected.";
}

TEST(StretchingFunctionTest, ConcentratesAtEndEvenWithInvertedInterval) {
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = true};
  StretchingFunction stretch(1.0, 0.0, 20, params); // Invertido
  auto nodes = stretch.generateNodes();

  // Expect nodes to be closer to each other at the end (0.0)
  double first_gap = std::abs(nodes[1] - nodes[0]); // near 1.0 (xBeg)
  double last_gap = std::abs(nodes[nodes.size() - 1] -
                             nodes[nodes.size() - 2]); // near 0.0 (xEnd)
  ASSERT_GT(first_gap, last_gap)
      << "Spacing near x=0.0 is not smaller, even with inverted range.";
}

TEST(StretchingFunctionTest, ConcentratesAtBeginningWithInvertedInterval) {
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = false};
  StretchingFunction stretch(1.0, 0.0, 20, params); // Invertido
  auto nodes = stretch.generateNodes();

  // Expect nodes to be closer to each other at the beginning (1.0)
  double first_gap = std::abs(nodes[1] - nodes[0]); // near 1.0 (xBeg)
  double last_gap = std::abs(nodes[nodes.size() - 1] -
                             nodes[nodes.size() - 2]); // near 0.0 (xEnd)
  ASSERT_LT(first_gap, last_gap)
      << "Spacing near x=1.0 is not smaller when stretching is reversed.";
}

TEST(StretchingFunctionTest, ThrowsOnInvalidRange) {
  StretchingParams params{};
  // Points cannot be the same
  EXPECT_THROW(StretchingFunction(1.0, 1.0, 10, params), std::invalid_argument);
  // Accept inverted domains
  EXPECT_NO_THROW(StretchingFunction(2.0, 1.0, 10, params));
}

TEST(StretchingFunctionTest, ThrowsOnInvalidNumPoints) {
  StretchingParams params{};
  EXPECT_THROW(StretchingFunction(0.0, 1.0, 1, params), std::invalid_argument);
  EXPECT_THROW(StretchingFunction(0.0, 1.0, 0, params), std::invalid_argument);
}

TEST(StretchingFunctionTest, NodesAreMonotonicIncreasing) {
  StretchingParams params{};
  StretchingFunction stretch(-2.0, 2.0, 50, params);
  std::vector<double> nodes = stretch.generateNodes();

  for (size_t i = 1; i < nodes.size(); ++i) {
    ASSERT_GT(nodes[i], nodes[i - 1])
        << "Node at index " << i << " is not greater than previous.";
  }
}

TEST(StretchingFunctionTest, UniformStretchingGivesLinearSpacing) {
  StretchingParams params{.paramP = 1.0, .paramQ = 2.0}; // P = 1 is uniform
  StretchingFunction stretch(0.0, 1.0, 5, params);
  auto nodes = stretch.generateNodes();

  double expected_dx = 1.0 / 4.0;
  for (int i = 0; i < 5; ++i) {
    ASSERT_NEAR(nodes[i], i * expected_dx, SOLVER_TOLERANCE);
  }
}

TEST(StretchingFunctionTest, DerivativesCorrectSize) {
  StretchingParams params{};
  StretchingFunction stretch(0.0, 1.0, 11, params);
  std::vector<double> derivatives = stretch.generateNodeDerivatives();

  ASSERT_EQ(derivatives.size(), 11)
      << "Derivatives vector size does not match numPoints.";
}

TEST(StretchingFunctionTest, DerivativesConstantForP1) {
  StretchingParams params{.paramP = 1.0, .paramQ = 2.0}; // p = 1 is uniform
  double xBeg = 0.0;
  double xEnd = 2.0;
  int numPoints = 5;
  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto derivatives = stretch.generateNodeDerivatives();

  double expected_dx_deta = xEnd - xBeg; // Constant derivative for p = 1
  for (size_t i = 0; i < derivatives.size(); ++i) {
    ASSERT_NEAR(derivatives[i], expected_dx_deta, SOLVER_TOLERANCE)
        << "Derivative at index " << i << " is not constant for p = 1.";
  }
}

TEST(StretchingFunctionTest, DerivativesArePositive) {
  StretchingParams params{};
  StretchingFunction stretch(-1.0, 1.0, 20, params);
  auto derivatives = stretch.generateNodeDerivatives();

  for (size_t i = 0; i < derivatives.size(); ++i) {
    ASSERT_GT(derivatives[i], 0.0)
        << "Derivative at index " << i << " is not positive.";
  }
}

TEST(StretchingFunctionTest, DerivativesMatchAnalyticalValues) {
  StretchingParams params{.paramP = 0.5, .paramQ = 2.0};
  double xBeg = 0.0;
  double xEnd = 1.0;
  int numPoints = 3; // Test at few points for simplicity
  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto derivatives = stretch.generateNodeDerivatives();

  // Pre-calculated analytical derivatives for eta = {0, 0.5, 1}
  double q = params.paramQ;
  double p = params.paramP;
  double dx = xEnd - xBeg;
  double tanhQ = std::tanh(q);
  std::vector<double> eta = {0.0, 0.5, 1.0};
  std::vector<double> expected_derivatives;

  for (double e : eta) {
    double sechSquared = 1.0 - std::pow(std::tanh(q * (1.0 - e)), 2);
    double ds_deta = p + ((1.0 - p) * q * sechSquared / tanhQ);
    expected_derivatives.push_back(dx * ds_deta);
  }

  for (size_t i = 0; i < derivatives.size(); ++i) {
    ASSERT_NEAR(derivatives[i], expected_derivatives[i], SOLVER_TOLERANCE)
        << "Derivative at index " << i << " does not match expected value.";
  }
}

TEST(StretchingFunctionTest, DerivativesAreNegativeWithInvertedDomain) {
  StretchingParams params{
      .paramP = 1.9, .paramQ = 2.0, .concentrateAtEnd = true};
  StretchingFunction stretch(1.0, 0.0, 30, params); // Inverted
  auto derivatives = stretch.generateNodeDerivatives();

  for (size_t i = 0; i < derivatives.size(); ++i) {
    ASSERT_LT(derivatives[i], 0.0) << "Derivative at index " << i
                                   << " is not negative with inverted domain.";
  }
}
