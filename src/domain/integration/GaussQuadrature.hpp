#pragma once

#include "IQuadratureRule.hpp"
#include <cmath>
#include <vector>

namespace bem::domain::integration {

class GaussQuadrature : public IQuadratureRule {
public:
  explicit GaussQuadrature(int order);

  [[nodiscard]] Complex integrate(
      const Panel &panel,
      const std::function<Complex(const Point2D &)> &integrand) const override;

  [[nodiscard]] int getOrder() const {
    return order_;
  }
  [[nodiscard]] const std::vector<QuadraturePoint> &getPoints() const {
    return points_;
  }

private:
  int order_;
  std::vector<QuadraturePoint> points_;

  void setupQuadraturePoints();
};

} // namespace bem::domain::integration
