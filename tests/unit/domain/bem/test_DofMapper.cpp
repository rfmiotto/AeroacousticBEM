#include "domain/bem/DofMapper.hpp"
#include <gtest/gtest.h>

using namespace bem::foundation::dof;
using namespace bem::types;

namespace {

TEST(DofMapperTest, EmptyMeshHasNoDofs) {
  std::vector<Element> elems;
  DofMapper mapper(elems);
  EXPECT_EQ(mapper.nGlobalDofs(), 0);
}

TEST(DofMapperTest, ConstantElementProducesMidpointAsGlobalDof) {
  Point2D p1{0.0, 0.0};
  Point2D p2{1.0, 0.0};
  std::vector<Element> elems = {Element(p1, p2, ElementType::CONSTANT)};
  DofMapper mapper(elems);

  ASSERT_EQ(mapper.nGlobalDofs(), 1);
  const auto &gp = mapper.globalPoints();
  EXPECT_NEAR(gp[0].x, 0.5, 1e-12);
  EXPECT_NEAR(gp[0].y, 0.0, 1e-12);

  const auto &ids = mapper.getElementGlobalIds(0);
  ASSERT_EQ(ids.size(), 1);
  EXPECT_EQ(ids[0], 0);
}

TEST(DofMapperTest, LinearTwoElementsSharesNode) {
  Point2D p1{0.0, 0.0};
  Point2D p2{1.0, 0.0};
  Point2D p3{2.0, 0.0};
  std::vector<Element> elems = {Element(p1, p2, ElementType::LINEAR),
                                Element(p2, p3, ElementType::LINEAR)};
  DofMapper mapper(elems);

  ASSERT_EQ(mapper.nGlobalDofs(), 3);
  const auto &gp = mapper.globalPoints();
  EXPECT_EQ(gp[0], p1);
  EXPECT_EQ(gp[1], p2);
  EXPECT_EQ(gp[2], p3);

  const auto &ids0 = mapper.getElementGlobalIds(0);
  ASSERT_EQ(ids0.size(), 2);
  EXPECT_EQ(ids0[0], 0);
  EXPECT_EQ(ids0[1], 1);

  const auto &ids1 = mapper.getElementGlobalIds(1);
  ASSERT_EQ(ids1.size(), 2);
  EXPECT_EQ(ids1[0], 1); // shared node
  EXPECT_EQ(ids1[1], 2);

  // Reverse mapping
  const auto &owners = mapper.getGlobalOwners(1);
  ASSERT_EQ(owners.size(), 2); // shared between both elements
  EXPECT_EQ(owners[0].first, 0);
  EXPECT_EQ(owners[1].first, 1);
}

TEST(DofMapperTest, QuadraticElementAddsMidpoint) {
  Point2D p1{0.0, 0.0};
  Point2D p2{2.0, 0.0};
  std::vector<Element> elems = {Element(p1, p2, ElementType::QUADRATIC)};
  DofMapper mapper(elems);

  ASSERT_EQ(mapper.nGlobalDofs(), 3);
  const auto &gp = mapper.globalPoints();
  EXPECT_EQ(gp[0], p1);
  EXPECT_EQ(gp[2], p2);
  EXPECT_NEAR(gp[1].x, 1.0, 1e-12); // midpoint
  EXPECT_NEAR(gp[1].y, 0.0, 1e-12);

  const auto &ids = mapper.getElementGlobalIds(0);
  ASSERT_EQ(ids.size(), 3);
  EXPECT_EQ(ids[0], 0);
  EXPECT_EQ(ids[1], 1);
  EXPECT_EQ(ids[2], 2);
}

TEST(DofMapperTest, ClosureMergesFirstAndLastNode) {
  Point2D p1{0.0, 0.0};
  Point2D p2{1.0, 0.0};
  Point2D p3{0.0, 0.0};
  std::vector<Element> elems = {Element(p1, p2, ElementType::LINEAR),
                                Element(p2, p3, ElementType::LINEAR)};
  DofMapper mapper(elems);

  // Closure merges node p3 with p1, so only two unique points remain
  ASSERT_EQ(mapper.nGlobalDofs(), 2);
  const auto &gp = mapper.globalPoints();
  EXPECT_EQ(gp[0], p1);
  EXPECT_EQ(gp[1], p2);

  // Check mapping
  const auto &ids0 = mapper.getElementGlobalIds(0);
  const auto &ids1 = mapper.getElementGlobalIds(1);
  EXPECT_EQ(ids0[0], 0);
  EXPECT_EQ(ids0[1], 1);
  EXPECT_EQ(ids1[0], 1);
  EXPECT_EQ(ids1[1], 0); // closed back to p1
}

TEST(DofMapperTest, ThrowsOnInvalidAccess) {
  Point2D p1{0.0, 0.0};
  Point2D p2{1.0, 0.0};
  std::vector<Element> elems = {Element(p1, p2, ElementType::LINEAR)};
  DofMapper mapper(elems);

  // Wrap calls in lambdas and assign the return value to a local variable
  // so [[nodiscard]] is satisfied at compile-time while still allowing the
  // exception to be thrown at runtime.
  EXPECT_THROW(
      [&]() {
        [[maybe_unused]] const auto &unused = mapper.getElementGlobalIds(10);
      }(),
      std::out_of_range);

  EXPECT_THROW(
      [&]() {
        [[maybe_unused]] const auto &unused = mapper.getGlobalOwners(10);
      }(),
      std::out_of_range);

  EXPECT_THROW(
      [&]() {
        [[maybe_unused]] auto unused = mapper.elementLocalToGlobal(0, 10);
      }(),
      std::out_of_range);

  EXPECT_THROW(
      [&]() {
        [[maybe_unused]] auto unused = mapper.elementLocalToGlobal(10, 0);
      }(),
      std::out_of_range);
}

} // namespace
