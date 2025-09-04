#include "domain/integration/AnalyticalConstantSelfTerm.hpp"
#include "domain/integration/GaussQuadrature.hpp"
#include "domain/integration/IntegrationStrategy.hpp"
#include "domain/integration/TellesQuadrature.hpp"
#include "gtest/gtest.h"

using namespace bem::domain::integration;
using namespace bem::types;

namespace {

// ---- isNearSingular tests ----

TEST(IsNearSingularTest, PointOnSegment) {
  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR};
  Point2D pt{0.5, 0.0}; // exactly on segment
  EXPECT_TRUE(isNearSingular(elem, pt, 1e-12));
}

TEST(IsNearSingularTest, PointCloseToSegment) {
  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR};
  Point2D pt{0.5, 1e-5}; // very close above
  EXPECT_TRUE(isNearSingular(elem, pt, 1e-4));
}

TEST(IsNearSingularTest, PointFarFromSegment) {
  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR};
  Point2D pt{0.5, 1.0}; // clearly away
  EXPECT_FALSE(isNearSingular(elem, pt, 1e-4));
}

// ---- selectQuadratureRule tests ----

TEST(SelectQuadratureRuleTest, AnalyticalConstantElementSelfTerm) {
  IntegrationParameters ip{};
  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::CONSTANT};
  Point2D colloc{0.5, 0.0}; // on element
  double k = 2.0;

  auto rule = selectQuadratureRule(ip, elem, colloc, k);

  ASSERT_NE(rule, nullptr);
  EXPECT_NE(dynamic_cast<AnalyticalConstantQuadrature *>(rule.get()), nullptr);
}

TEST(SelectQuadratureRuleTest, NearSingularUsesTelles) {
  IntegrationParameters ip{};
  ip.singularity_distance = 0.1;

  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR};
  Point2D colloc{0.5, 1e-2}; // near but not on element
  double k = 2.0;

  auto rule = selectQuadratureRule(ip, elem, colloc, k);

  ASSERT_NE(rule, nullptr);
  EXPECT_NE(dynamic_cast<TellesQuadrature *>(rule.get()), nullptr);
}

TEST(SelectQuadratureRuleTest, DefaultGaussQuadrature) {
  IntegrationParameters ip{};
  ip.singularity_distance = 1e-3;

  Element elem{Point2D{0.0, 0.0}, Point2D{1.0, 0.0}, ElementType::LINEAR};
  Point2D colloc{0.5, 1.0}; // far away
  double k = 2.0;

  auto rule = selectQuadratureRule(ip, elem, colloc, k);

  ASSERT_NE(rule, nullptr);
  EXPECT_NE(dynamic_cast<GaussQuadrature *>(rule.get()), nullptr);
}

} // namespace
