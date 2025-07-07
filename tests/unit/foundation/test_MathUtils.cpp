#include <cmath>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numbers>

#include "foundation/types/MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"
#include "foundation/utils/MathUtils.hpp"

using namespace bem::foundation;
using namespace bem::types;
using namespace bem::foundation::utils;
using ::testing::DoubleNear;

class MathUtilsTest : public ::testing::Test {
protected:
  double tolerance_;
  Point2D point1_, point2_, point3_;

  void SetUp() override {
    // Common test data
    tolerance_ = 1e-12;
    point1_ = Point2D{0.0, 0.0};
    point2_ = Point2D{3.0, 4.0};
    point3_ = Point2D{1.0, 1.0};
  }
};

// Test Constants
TEST(ConstantsTest, MathematicalConstants) {
  EXPECT_THAT(Constants::PI, DoubleNear(std::numbers::pi, 1e-15));
  EXPECT_THAT(Constants::E, DoubleNear(std::numbers::e, 1e-15));
  EXPECT_THAT(Constants::EULER_MASCHERONI,
              DoubleNear(0.5772156649015329, 1e-15));
  EXPECT_THAT(Constants::SQRT_PI,
              DoubleNear(std::sqrt(std::numbers::pi), 1e-15));
  EXPECT_THAT(Constants::INV_4PI,
              DoubleNear(1.0 / (4.0 * std::numbers::pi), 1e-15));
}

TEST(ConstantsTest, NumericalConstants) {
  EXPECT_GT(Constants::MACHINE_EPSILON, 0.0);
  EXPECT_LT(Constants::MACHINE_EPSILON, 1e-10);
  EXPECT_GT(Constants::INTEGRATION_TOLERANCE, Constants::MACHINE_EPSILON);
  EXPECT_GT(Constants::SOLVER_TOLERANCE, Constants::MACHINE_EPSILON);
}

// Test basic mathematical utilities
TEST_F(MathUtilsTest, DistanceCalculation) {
  double expected_distance = 5.0; // 3-4-5 triangle
  EXPECT_THAT(MathUtils::distance(point1_, point2_),
              DoubleNear(expected_distance, tolerance_));

  // Distance should be symmetric
  EXPECT_THAT(MathUtils::distance(point2_, point1_),
              DoubleNear(expected_distance, tolerance_));

  // Distance to self should be zero
  EXPECT_THAT(MathUtils::distance(point1_, point1_),
              DoubleNear(0.0, tolerance_));
}

TEST_F(MathUtilsTest, DistanceSquaredCalculation) {
  double expected_distance_sq = 25.0; // 5^2
  EXPECT_THAT(MathUtils::distanceSquared(point1_, point2_),
              DoubleNear(expected_distance_sq, tolerance_));

  // Should be more efficient than computing distance and squaring
  double dist = MathUtils::distance(point1_, point2_);
  EXPECT_THAT(MathUtils::distanceSquared(point1_, point2_),
              DoubleNear(dist * dist, tolerance_));
}

TEST_F(MathUtilsTest, AngleCalculation) {
  Point2D origin{0.0, 0.0};
  Point2D positive_x{1.0, 0.0};
  Point2D positive_y{0.0, 1.0};
  Point2D negative_x{-1.0, 0.0};

  EXPECT_NEAR(MathUtils::angle(origin, positive_x), 0.0, tolerance_);
  EXPECT_NEAR(MathUtils::angle(origin, positive_y), Constants::PI / 2.0,
              tolerance_);
  EXPECT_NEAR(MathUtils::angle(origin, negative_x), Constants::PI, tolerance_);
}

TEST_F(MathUtilsTest, IsZeroFunction) {
  EXPECT_TRUE(MathUtils::isZero(0.0));
  EXPECT_TRUE(MathUtils::isZero(1e-16));
  EXPECT_TRUE(MathUtils::isZero(-1e-16));
  EXPECT_FALSE(MathUtils::isZero(1e-10));
  EXPECT_FALSE(MathUtils::isZero(-1e-10));

  // Custom tolerance
  EXPECT_TRUE(MathUtils::isZero(1e-5, 1e-4));
  EXPECT_FALSE(MathUtils::isZero(1e-3, 1e-4));
}

TEST_F(MathUtilsTest, IsEqualFunction) {
  EXPECT_TRUE(MathUtils::isEqual(1.0, 1.0));
  EXPECT_TRUE(MathUtils::isEqual(1.0, 1.0 + 1e-16));
  EXPECT_FALSE(MathUtils::isEqual(1.0, 1.1));

  // Custom tolerance
  EXPECT_TRUE(MathUtils::isEqual(1.0, 1.01, 0.02));
  EXPECT_FALSE(MathUtils::isEqual(1.0, 1.01, 0.005));
}

// Test complex number utilities
TEST_F(MathUtilsTest, ComplexMagnitudeAndPhase) {
  Complex z1{3.0, 4.0};
  Complex z2{-1.0, 1.0};
  Complex z3{0.0, 0.0};

  EXPECT_THAT(MathUtils::magnitude(z1), DoubleNear(5.0, tolerance_));
  EXPECT_THAT(MathUtils::magnitude(z2), DoubleNear(std::sqrt(2.0), tolerance_));
  EXPECT_THAT(MathUtils::magnitude(z3), DoubleNear(0.0, tolerance_));

  EXPECT_THAT(MathUtils::phase(z1),
              DoubleNear(std::atan2(4.0, 3.0), tolerance_));
  EXPECT_THAT(MathUtils::phase(z2),
              DoubleNear(3.0 * Constants::PI / 4.0, tolerance_));
}
