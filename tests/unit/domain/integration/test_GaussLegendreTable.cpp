#include "domain/integration/GaussLegendreTable.hpp"
#include "foundation/exceptions/BEMException.hpp"
#include "foundation/utils/Constants.hpp"
#include <gtest/gtest.h>

using namespace bem::domain::integration;
using namespace bem::foundation::exceptions;
using namespace bem::foundation::utils::Constants;

TEST(GaussLegendreTableTest, CorrectNumberOfPoints) {
  std::vector<int> orders = {1, 2, 3, 4, 5, 6, 8};
  for (int order : orders) {
    auto pts = GaussLegendreTable::getPoints(order);
    EXPECT_EQ(static_cast<int>(pts.size()), order)
        << "Unexpected number of points for order " << order;
  }
}

TEST(GaussLegendreTableTest, WeightsSumToTwo) {
  std::vector<int> orders = {1, 2, 3, 4, 5, 6, 8};
  for (int order : orders) {
    auto pts = GaussLegendreTable::getPoints(order);
    double sum_w = 0.0;
    for (auto &qp : pts) {
      sum_w += qp.weight;
    }
    EXPECT_NEAR(sum_w, 2.0, INTEGRATION_TOLERANCE)
        << "Weight sum incorrect for order " << order;
  }
}

TEST(GaussLegendreTableTest, SymmetryOfNodesAndWeights) {
  std::vector<int> orders = {2, 3, 4,
                             5, 6, 8}; // skip order=1 (only one point at 0)
  for (int order : orders) {
    auto pts = GaussLegendreTable::getPoints(order);

    for (int i = 0; i < order / 2; ++i) {
      double xi_left = pts[i].xi;
      double xi_right = pts[order - 1 - i].xi;
      double w_left = pts[i].weight;
      double w_right = pts[order - 1 - i].weight;

      EXPECT_NEAR(xi_left, -xi_right, INTEGRATION_TOLERANCE)
          << "Symmetry mismatch for xi at order " << order;
      EXPECT_NEAR(w_left, w_right, INTEGRATION_TOLERANCE)
          << "Symmetry mismatch for weight at order " << order;
    }
  }
}

TEST(GaussLegendreTableTest, UnsupportedOrderThrows) {
  EXPECT_THROW(GaussLegendreTable::getPoints(7), BEMIntegrationException);
  EXPECT_THROW(GaussLegendreTable::getPoints(0), BEMIntegrationException);
  EXPECT_THROW(GaussLegendreTable::getPoints(10), BEMIntegrationException);
}
