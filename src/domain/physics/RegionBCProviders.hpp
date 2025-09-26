#pragma once

#include "BoundaryConditionProvider.hpp"
#include "types/BoundaryValueTypes.hpp"

namespace bem::domain::bem {

using types::BoundaryCondition;
using types::Complex;
using types::Real;

/**
 * @brief Upstream region: Dirichlet = fixed value (e.g., 0)
 */
class UpstreamBCProvider : public IBoundaryConditionProvider {
public:
  [[nodiscard]] BoundaryCondition evaluate(Real /*x*/) const override {
    return BoundaryCondition{types::BCType::DIRICHLET, Complex(0.0, 0.0)};
  }

  [[nodiscard]] types::BCType type() const override {
    return types::BCType::DIRICHLET;
  }
};

/**
 * @brief Airfoil region: Dirichlet = function of x
 */
class AirfoilBCProvider : public IBoundaryConditionProvider {
public:
  AirfoilBCProvider(double w, double beta, double k1)
      : w_(w), beta_(beta), k1_(k1) {
  }

  [[nodiscard]] BoundaryCondition evaluate(double x) const override {
    using namespace std::complex_literals;
    std::complex<double> val =
        (w_ / (2.0 * beta_)) * std::exp(-1i * k1_ * x / (2.0 * beta_ * beta_));

    return BoundaryCondition{types::BCType::NEUMANN, val};
  }

  [[nodiscard]] types::BCType type() const override {
    return types::BCType::NEUMANN;
  }

private:
  double w_, beta_, k1_;
};

/**
 * @brief Downstream region: Kutta condition (requires neighbors)
 */
class DownstreamBCProvider : public IBoundaryConditionProvider {
public:
  DownstreamBCProvider(double w, double beta, double k1)
      : w_(w), beta_(beta), k1_(k1) {
  }

  [[nodiscard]] BoundaryCondition evaluate(double /*x*/) const override {
    // TODO: Proper implementation with x-derivatives
    double dummy = w_ + beta_ + k1_;
    return BoundaryCondition{types::BCType::KUTTA,
                             std::complex<double>(dummy, 0.0)};
  }

  [[nodiscard]] types::BCType type() const override {
    return types::BCType::KUTTA;
  }

private:
  double w_, beta_, k1_;
};

} // namespace bem::domain::bem
