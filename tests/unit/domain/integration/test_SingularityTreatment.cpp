#include "domain/integration/SingularityTreatment.hpp"
#include "domain/integration/TellesQuadrature.hpp"
#include "foundation/utils/Constants.hpp"
#include "integration/GreenFunctionIntegrator.hpp"
#include <gtest/gtest.h>

using namespace bem::domain::integration;
using namespace bem::types;
using namespace bem::foundation::utils::Constants;

namespace {
Element createUnitPanel() {
  return {Point2D{0.0, 0.0}, Point2D{1.0, 0.0}};
}
} // namespace

TEST(ConstantSingularityTreatmentTest, GreenFunctionAnalytical) {
  ConstantSingularityTreatment treatment;
  Element element = createUnitPanel();
  double k = 1.0;

  Complex result = treatment.treatGreen(element, element.midpoint(), k);

  double expected_real =
      (1.0 - EULER_MASCHERONI - std::log(k * 0.25)) * INV_2PI;
  double expected_imag = 0.25;

  EXPECT_NEAR(result.real(), expected_real, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result.imag(), expected_imag, INTEGRATION_TOLERANCE);
}

TEST(ConstantSingularityTreatmentTest, NormalDerivativeAnalytical) {
  ConstantSingularityTreatment treatment;
  Element element = createUnitPanel();
  Point2D midpoint = element.midpoint();
  Vector2D normal{0.0, 1.0}; // Arbitrary

  Complex result =
      treatment.treatNormalDerivative(element, midpoint, normal, 1.0);
  EXPECT_NEAR(result.real(), 0.5, INTEGRATION_TOLERANCE);
  EXPECT_NEAR(result.imag(), 0.0, INTEGRATION_TOLERANCE);
}
