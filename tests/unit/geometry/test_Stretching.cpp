#include "domain/geometry/Stretching.hpp"
#include <gtest/gtest.h>

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

  ASSERT_NEAR(nodes.front(), xBeg, 1e-12);
  ASSERT_NEAR(nodes.back(), xEnd, 1e-12);
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
    ASSERT_NEAR(nodes[i], i * expected_dx, 1e-12);
  }
}
