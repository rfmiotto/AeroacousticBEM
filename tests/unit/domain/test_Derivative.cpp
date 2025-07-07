#include "domain/derivative/Derivative.hpp"
#include "foundation/utils/Constants.hpp"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

using namespace bem::foundation::utils::Constants;

template <typename T> class CompactScheme10Fixture : public ::testing::Test {
public:
  using Vec = std::vector<T>;

  static constexpr std::size_t N = 100;
  static constexpr double A = 0.0;
  static constexpr double B = 2.0 * PI;
  const double h = (B - A) / (N - 1);

  std::vector<double> x;
  Vec f, expected;

  // f(x) = sin(x), f'(x) = cos(x)
  void SetUp() override {
    x = linspace(A, B, N);
    f.resize(N);
    expected.resize(N);
    for (std::size_t i = 0; i < N; ++i) {
      f[i] = std::sin(x[i]);
      expected[i] = std::cos(x[i]);
    }
  }

private:
  // Helper function to generate equally spaced points
  static std::vector<double> linspace(double start, double end, std::size_t n) {
    std::vector<double> result(n);
    double step = (end - start) / static_cast<double>(n - 1);
    for (std::size_t i = 0; i < n; ++i) {
      result[i] = start + static_cast<double>(i) * step;
    }
    return result;
  }
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(CompactScheme10Fixture, MyTypes);

TYPED_TEST(CompactScheme10Fixture, SineFunctionDerivative) {
  using T = TypeParam;

  CompactScheme10<T> scheme(this->N);
  auto numerical = scheme.computeDerivative(this->f, this->h);

  // Save to CSV
  // std::ofstream out(std::is_same_v<T, double>
  //                       ? "derivative_results_real.csv"
  //                       : "derivative_results_complex.csv");
  //
  // out << "x,numerical,expected\n";
  // for (std::size_t i = 0; i < this->N; ++i) {
  //   if constexpr (std::is_same_v<T, std::complex<double>>) {
  //     out << this->x[i] << "," << numerical[i].real() << ","
  //         << this->expected[i].real() << "\n";
  //   } else {
  //     out << this->x[i] << "," << numerical[i] << "," << this->expected[i]
  //         << "\n";
  //   }
  // }
  // out.close();

  // Allowable tolerance
  const double tol = 5e-3;

  // Ignore boundary values (less accurate)
  for (std::size_t i = 5; i < this->N - 5; ++i) {
    if constexpr (std::is_same_v<T, std::complex<double>>) {
      EXPECT_NEAR(numerical[i].real(), this->expected[i].real(), tol);
    } else {
      EXPECT_NEAR(numerical[i], this->expected[i], tol);
    }
  }
}

TYPED_TEST(CompactScheme10Fixture, MultipleCalls) {
  using T = TypeParam;

  CompactScheme10<T> scheme(this->N);
  auto firstDerivative = scheme.computeDerivative(this->f, this->h);
  auto secondDerivative = scheme.computeDerivative(firstDerivative, this->h);

  // Save to CSV
  // std::ofstream out(std::is_same_v<T, double>
  //                       ? "derivative_results_real.csv"
  //                       : "derivative_results_complex.csv");
  //
  // out << "x,numerical,expected\n";
  // for (std::size_t i = 0; i < this->N; ++i) {
  //   if constexpr (std::is_same_v<T, std::complex<double>>) {
  //     out << this->x[i] << "," << secondDerivative[i].real() << ","
  //         << -this->expected[i].real() << "\n";
  //   } else {
  //     out << this->x[i] << "," << secondDerivative[i] << "," << -this->f[i]
  //         << "\n";
  //   }
  // }
  // out.close();

  // Allowable tolerance
  const double tol = 5e-3;

  // Ignore boundary values (less accurate)
  for (std::size_t i = 5; i < this->N - 5; ++i) {
    if constexpr (std::is_same_v<T, std::complex<double>>) {
      EXPECT_NEAR(secondDerivative[i].real(), -this->f[i].real(), tol);
    } else {
      EXPECT_NEAR(secondDerivative[i], -this->f[i], tol);
    }
  }
}
