#include "domain/integration/GreenFunctionIntegrator.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace bem::domain::integration;
using namespace bem::types;
using namespace testing;

class MockSingularityTreatment : public SingularityTreatment {
public:
  MOCK_METHOD(Complex,
              treatGreen,
              (const Element &, const Point2D &, double),
              (const, override));
  MOCK_METHOD(Complex,
              treatNormalDerivative,
              (const Element &, const Point2D &, const Vector2D &, double),
              (const, override));
};

// Case 1: singularity + constant element → use singularity_handler
TEST(GreenFunctionIntegratorTest, ConstantSingularElementUsesHandler) {
  IntegrationParameters params;
  params.order = 3;
  params.use_singularity_treatment = true;

  Element elem({0, 0}, {1, 0});
  elem.type = ElementType::CONSTANT;

  Point2D x = elem.midpoint(); // is over the element → singularidade

  auto mock_handler = std::make_shared<MockSingularityTreatment>();
  EXPECT_CALL(*mock_handler, treatGreen(_, _, _))
      .WillOnce(Return(Complex{42.0, 0.0}));

  auto integrator = makeGreenIntegrator(params, mock_handler);
  Complex result = integrator(elem, x, 1.0);

  EXPECT_EQ(result, Complex(42.0, 0.0));
}

// Case 2: singularity + non-constant + use_singularity_treatment = true →
// use TellesQuadrature
TEST(GreenFunctionIntegratorTest, SingularNonConstantUsesTelles) {
  IntegrationParameters params;
  params.order = 2;
  params.use_singularity_treatment = true;

  Element elem({0, 0}, {1, 0});
  elem.type = ElementType::LINEAR; // não constante

  Point2D x = elem.midpoint(); // force contains=true

  auto integrator = makeGreenIntegrator(params);
  Complex result = integrator(elem, x, 0.5);

  // Only checks for finite value
  EXPECT_TRUE(std::isfinite(result.real()));
  EXPECT_TRUE(std::isfinite(result.imag()));
}

// Case 3: singularity + non-constant + use_singularity_treatment = false →
// use GaussQuadrature
TEST(GreenFunctionIntegratorTest, SingularNonConstantNoTreatmentUsesGauss) {
  IntegrationParameters params;
  params.order = 2;
  params.use_singularity_treatment = false;

  Element elem({0, 0}, {1, 0});
  elem.type = ElementType::LINEAR;

  Point2D x = elem.midpoint();

  auto integrator = makeGreenIntegrator(params);
  Complex result = integrator(elem, x, 0.5);

  EXPECT_TRUE(std::isfinite(result.real()));
  EXPECT_TRUE(std::isfinite(result.imag()));
}

// Case 4: non-singular → GaussQuadrature
TEST(GreenFunctionIntegratorTest, NonSingularUsesGauss) {
  IntegrationParameters params;
  params.order = 2;
  params.use_singularity_treatment = true; // doesn't matter

  Element elem({0, 0}, {1, 0});
  elem.type = ElementType::LINEAR;

  Point2D x = {0.0, 5.0}; // far from element → contains=false

  auto integrator = makeGreenIntegrator(params);
  Complex result = integrator(elem, x, 0.5);

  EXPECT_TRUE(std::isfinite(result.real()));
  EXPECT_TRUE(std::isfinite(result.imag()));
}
