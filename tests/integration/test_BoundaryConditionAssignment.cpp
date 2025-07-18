#include "domain/geometry/FlatPlate.hpp"
#include "domain/physics/BoundaryConditionManager.hpp"
#include <gtest/gtest.h>

using namespace bem::boundary;
using namespace bem::domain::geometry;
using namespace bem::types;

namespace {

bem::types::Complex incidentPotential(const bem::types::Point2D &pt) {
  // Mock simples: valor arbitrário determinístico com base na posição
  return {pt.x + pt.y, pt.x - pt.y};
}

} // namespace

TEST(BoundaryIntegrationTest, DirichletConditionsAreCorrectlyAssigned) {
  FlatPlate2D plate(10, 1, 0.1, 3.0);

  BoundaryConditionManager bc_manager;

  for (ElementId id = 0; id < plate.getElements().size(); ++id) {
    const auto &element = plate.getElements()[id];

    if (element.region == ElementRegion::AIRFOIL) {
      Complex phi = incidentPotential(element.midpoint());
      bc_manager.setDirichlet(id, phi);
    }
  }

  for (ElementId id = 0; id < plate.getElements().size(); ++id) {
    const auto &element = plate.getElements()[id];
    if (element.region == ElementRegion::AIRFOIL) {
      ASSERT_TRUE(bc_manager.has(id));
      EXPECT_EQ(bc_manager.getType(id), BoundaryConditionType::DIRICHLET);

      Complex expected = incidentPotential(element.midpoint());
      EXPECT_EQ(bc_manager.getValue(id).potential, expected);
    } else {
      EXPECT_FALSE(bc_manager.has(id));
    }
  }
}
