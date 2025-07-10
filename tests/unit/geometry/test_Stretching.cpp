#include "domain/geometry/Stretching.hpp"
#include "utils/Constants.hpp"
#include <gtest/gtest.h>

using namespace bem::foundation::utils::Constants;

TEST(StretchingFunctionTest, GeneratesCorrectNumberOfNodes) {
  StretchingParams params{.paramP = 0.5, .paramQ = 3.0};
  StretchingFunction stretch(0.0, 1.0, 11, params);
  std::vector<double> nodes = stretch.generateNodes();

  ASSERT_EQ(nodes.size(), 11);
}

TEST(StretchingFunctionTest, NodesRespectDomainBounds) {
  StretchingParams params{.paramP = 0.7, .paramQ = 2.5};
  double xBeg = 0.0;
  double xEnd = 1.0;
  int numPoints = 11;

  StretchingFunction stretch(xBeg, xEnd, numPoints, params);
  auto nodes = stretch.generateNodes();

  ASSERT_NEAR(nodes.front(), xBeg, SOLVER_TOLERANCE);
  ASSERT_NEAR(nodes.back(), xEnd, SOLVER_TOLERANCE);
}

TEST(StretchingFunctionTest, ThrowsOnInvalidRange) {
  StretchingParams params{};
  EXPECT_THROW(StretchingFunction(1.0, 1.0, 10, params), std::invalid_argument);
  EXPECT_THROW(StretchingFunction(2.0, 1.0, 10, params), std::invalid_argument);
}

TEST(StretchingFunctionTest, ThrowsOnInvalidNumPoints) {
  StretchingParams params{};
  EXPECT_THROW(StretchingFunction(0.0, 1.0, 1, params), std::invalid_argument);
  EXPECT_THROW(StretchingFunction(0.0, 1.0, 0, params), std::invalid_argument);
}

TEST(StretchingFunctionTest, NodesAreMonotonicIncreasing) {
  StretchingParams params{.paramP = 0.3, .paramQ = 4.0};
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
  StretchingParams params{.paramP = 0.5, .paramQ = 3.0};
  StretchingFunction stretch(0.0, 1.0, 11, params);
  std::vector<double> derivatives = stretch.generateNodeDerivatives();

  ASSERT_EQ(derivatives.size(), 11)
      << "Derivatives vector size does not match numPoints.";
}

TEST(StretchingFunctionTest, DerivativesConstantForP1) {
  StretchingParams params{.paramP = 1.0,
                          .paramQ = 2.0}; // p = 1 gives linear spacing
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
  StretchingParams params{.paramP = 0.3, .paramQ = 4.0};
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
