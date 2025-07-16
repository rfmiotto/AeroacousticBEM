#include "domain/physics/BoundaryConditionManager.hpp"
#include "foundation/exceptions/BEMException.hpp"
#include <gtest/gtest.h>

using namespace bem::boundary;
using namespace bem::types;
using namespace bem::foundation::exceptions;

class BoundaryConditionManagerTest : public ::testing::Test {
protected:
  BoundaryConditionManager manager_;
  const ElementId id_ = 42;
};

TEST_F(BoundaryConditionManagerTest, SetDirichlet_StoresCorrectly) {
  Complex phi(1.0, 2.0);
  manager_.setDirichlet(id_, phi);

  EXPECT_TRUE(manager_.has(id_));
  EXPECT_EQ(manager_.getType(id_), BoundaryConditionType::DIRICHLET);
  EXPECT_EQ(manager_.getValue(id_).potential, phi);
  EXPECT_EQ(manager_.getValue(id_).normal_derivative, Complex(0.0));
}

TEST_F(BoundaryConditionManagerTest, SetNeumann_StoresCorrectly) {
  Complex dphi(0.5, -0.3);
  manager_.setNeumann(id_, dphi);

  EXPECT_TRUE(manager_.has(id_));
  EXPECT_EQ(manager_.getType(id_), BoundaryConditionType::NEUMANN);
  EXPECT_EQ(manager_.getValue(id_).normal_derivative, dphi);
  EXPECT_EQ(manager_.getValue(id_).potential, Complex(0.0));
}

TEST_F(BoundaryConditionManagerTest, SetRobin_StoresCorrectly) {
  Complex phi(0.1, 0.2);
  Complex dphi(0.3, 0.4);
  manager_.setRobin(id_, phi, dphi);

  EXPECT_TRUE(manager_.has(id_));
  EXPECT_EQ(manager_.getType(id_), BoundaryConditionType::ROBIN);
  EXPECT_EQ(manager_.getValue(id_).potential, phi);
  EXPECT_EQ(manager_.getValue(id_).normal_derivative, dphi);
}

TEST_F(BoundaryConditionManagerTest, Clear_RemovesEntry) {
  manager_.setDirichlet(id_, Complex(1.0));
  manager_.clear(id_);

  EXPECT_FALSE(manager_.has(id_));
}

TEST_F(BoundaryConditionManagerTest, Has_ReturnsTrueIfExists) {
  manager_.setNeumann(id_, Complex(2.0));
  EXPECT_TRUE(manager_.has(id_));
}

TEST_F(BoundaryConditionManagerTest, Has_ReturnsFalseIfNotExists) {
  EXPECT_FALSE(manager_.has(id_));
}

TEST_F(BoundaryConditionManagerTest, GetType_ThrowsIfNotFound) {

  EXPECT_THROW(
      { static_cast<void>(manager_.getType(id_)); },
      BoundaryConditionException);
}

TEST_F(BoundaryConditionManagerTest, GetValue_ThrowsIfNotFound) {
  EXPECT_THROW(
      { static_cast<void>(manager_.getValue(id_)); },
      BoundaryConditionException);
}
