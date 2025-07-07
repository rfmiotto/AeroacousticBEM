
#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

/**
 * @class CompactScheme10
 * @brief Computes the first-order derivative using a 10th-order compact finite
 * difference scheme.
 *
 * This class implements a high-order compact finite difference method to
 * compute the first-order derivative of a vector field. The scheme is based on
 * the 10th-order compact formulation described by Britton Olson in his Ph.D.
 * dissertation (Stanford University, 2012).
 *
 * @note See: Olson, B.J. (2012). *"Large-eddy simulation of multi-material
 * mixing and over-expanded nozzle flow"*, Ph.D. thesis, Stanford University.
 */
template <typename T> class CompactScheme10 {
public:
  using Vec = std::vector<T>;

  explicit CompactScheme10(std::size_t n) : n_(n) {
    if (n_ < 7) {
      throw std::invalid_argument(
          "Vector size must be at least 7 for 10th-order scheme.");
    }
    initializeDiagonals();
  }

  Vec computeDerivative(const Vec &vec, double h = 1.0) const {
    if (vec.size() != n_) {
      throw std::invalid_argument(
          "Input vector size does not match initialized size.");
    }

    Vec rhs(n_);

    rhs[0] = (-2.5 * vec[0] + 2.0 * vec[1] + 0.5 * vec[2]) / h;
    rhs[1] = (vec[2] - vec[0]) * A2 / h;
    rhs[2] = (vec[3] - vec[1]) * A3 / h + (vec[4] - vec[0]) * B3 / h;
    rhs[3] = (vec[4] - vec[2]) * A4 / h + (vec[5] - vec[1]) * B4 / h +
             (vec[6] - vec[0]) * C4 / h;

    std::vector<double> a_diag = A_;
    std::vector<double> b_diag = B_;
    std::vector<double> c_diag = C_;
    std::vector<double> d_diag = D_;
    std::vector<double> e_diag = E_;
    std::vector<double> f_diag = F_;

    for (std::size_t i = 4; i < n_ - 4; ++i) {
      rhs[i] = (vec[i + 1] - vec[i - 1]) * A / h +
               (vec[i + 2] - vec[i - 2]) * B / h +
               (vec[i + 3] - vec[i - 3]) * C / h;
    }

    rhs[n_ - 1] =
        (2.5 * vec[n_ - 1] - 2.0 * vec[n_ - 2] - 0.5 * vec[n_ - 3]) / h;
    rhs[n_ - 2] = (vec[n_ - 1] - vec[n_ - 3]) * A2 / h;
    rhs[n_ - 3] = (vec[n_ - 2] - vec[n_ - 4]) * A3 / h +
                  (vec[n_ - 1] - vec[n_ - 5]) * B3 / h;
    rhs[n_ - 4] = (vec[n_ - 3] - vec[n_ - 5]) * A4 / h +
                  (vec[n_ - 2] - vec[n_ - 6]) * B4 / h +
                  (vec[n_ - 1] - vec[n_ - 7]) * C4 / h;

    for (std::size_t i = 1; i < n_ - 1; ++i) {
      double m1 = a_diag[i - 1] / d_diag[i - 1];
      d_diag[i] -= m1 * c_diag[i - 1];
      c_diag[i] -= m1 * f_diag[i - 1];
      rhs[i] -= m1 * rhs[i - 1];

      double m2 = e_diag[i - 1] / d_diag[i - 1];
      a_diag[i] -= m2 * c_diag[i - 1];
      d_diag[i + 1] -= m2 * f_diag[i - 1];
      rhs[i + 1] -= m2 * rhs[i - 1];
    }

    double m = a_diag[n_ - 2] / d_diag[n_ - 2];
    d_diag[n_ - 1] -= m * c_diag[n_ - 2];

    Vec sol(n_);

    sol[n_ - 1] = (rhs[n_ - 1] - m * rhs[n_ - 2]) / d_diag[n_ - 1];
    sol[n_ - 2] = (rhs[n_ - 2] - c_diag[n_ - 2] * sol[n_ - 1]) / d_diag[n_ - 2];

    for (int i = static_cast<int>(n_) - 3; i >= 0; --i) {
      sol[i] = (rhs[i] - c_diag[i] * sol[i + 1] - f_diag[i] * sol[i + 2]) /
               d_diag[i];
    }

    return sol;
  }

private:
  std::size_t n_;
  std::vector<double> A_, B_, C_, D_, E_, F_;

  // Coeficientes do esquema
  static constexpr double ALPHA1 = 2.0;
  static constexpr double ALPHA2 = 0.25;
  static constexpr double A2 = 3.0 / 4.0;

  static constexpr double ALPHA3 = 4.7435 / 10.67175;
  static constexpr double BETA3 = 0.2964375 / 10.67175;
  static constexpr double A3 = 7.905 / 10.67175;
  static constexpr double B3 = 1.23515625 / 10.67175;

  static constexpr double ALPHA4 = 4.63271875 / 9.38146875;
  static constexpr double BETA4 = 0.451390625 / 9.38146875;
  static constexpr double A4 = 6.66984375 / 9.38146875;
  static constexpr double B4 = 1.53 / 9.38146875;
  static constexpr double C4 = 0.015 / 9.38146875;

  static constexpr double ALPHA = 0.5;
  static constexpr double BETA = 0.05;
  static constexpr double A = 17.0 / 24.0;
  static constexpr double B = 101.0 / 600.0;
  static constexpr double C = 0.01 / 6.0;

  void initializeDiagonals() {
    A_.assign(n_, ALPHA);
    B_.assign(n_, 0.0); // unused
    C_.assign(n_, ALPHA);
    D_.assign(n_, 1.0);
    E_.assign(n_, BETA);
    F_.assign(n_, BETA);

    // Ajustes para as fronteiras
    C_[0] = ALPHA1;
    C_[1] = ALPHA2;
    C_[2] = ALPHA3;
    C_[3] = ALPHA4;

    F_[0] = 0.0;
    F_[1] = 0.0;
    F_[2] = BETA3;
    F_[3] = BETA4;

    A_[0] = ALPHA2;
    A_[1] = ALPHA3;
    A_[2] = ALPHA4;

    E_[0] = BETA3;
    E_[1] = BETA4;

    E_[n_ - 1] = 0.0;
    E_[n_ - 2] = 0.0;
    E_[n_ - 3] = 0.0;
    E_[n_ - 4] = 0.0;
    E_[n_ - 5] = BETA3;
    E_[n_ - 6] = BETA4;

    A_[n_ - 1] = 0.0;
    A_[n_ - 2] = ALPHA1;
    A_[n_ - 3] = ALPHA2;
    A_[n_ - 4] = ALPHA3;
    A_[n_ - 5] = ALPHA4;

    C_[n_ - 1] = 0.0;
    C_[n_ - 2] = ALPHA2;
    C_[n_ - 3] = ALPHA3;
    C_[n_ - 4] = ALPHA4;

    F_[n_ - 1] = 0.0;
    F_[n_ - 2] = 0.0;
    F_[n_ - 3] = BETA3;
    F_[n_ - 4] = BETA4;
  }
};
