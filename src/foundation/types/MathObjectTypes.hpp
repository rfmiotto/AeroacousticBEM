#pragma once

#include <Eigen/Core>
#include <complex>
#include <vector>

namespace bem::types {

// Basic types
using Real = double;
using Complex = std::complex<Real>;

// Vector types
using RealVector = std::vector<Real>;
using ComplexVector = std::vector<Complex>;

// Matrix types (using Eigen for high-performance linear algebra)
using RealMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;
using RealVectorEigen = Eigen::VectorXd;
using ComplexVectorEigen = Eigen::VectorXcd;

// 2D geometric types
struct Point2D {
  Real x{0.0};
  Real y{0.0};

  Point2D() = default;
  Point2D(Real x_val, Real y_val) : x(x_val), y(y_val) {}

  // Arithmetic operations
  Point2D operator+(const Point2D &other) const {
    return {x + other.x, y + other.y};
  }

  Point2D operator-(const Point2D &other) const {
    return {x - other.x, y - other.y};
  }

  Point2D operator*(Real scalar) const { return {x * scalar, y * scalar}; }

  Point2D operator/(Real scalar) const { return {x / scalar, y / scalar}; }

  // Comparison operators
  bool operator==(const Point2D &other) const;
  bool operator!=(const Point2D &other) const;
};

// Vector2D for directions and normals
struct Vector2D {
  Real x{0.0};
  Real y{0.0};

  Vector2D() = default;
  Vector2D(Real x_val, Real y_val) : x(x_val), y(y_val) {}
  Vector2D(const Point2D &from, const Point2D &to)
      : x(to.x - from.x), y(to.y - from.y) {}

  // Vector operations
  [[nodiscard]] Real dot(const Vector2D &other) const {
    return (x * other.x) + (y * other.y);
  }

  [[nodiscard]] Real cross(const Vector2D &other) const {
    return (x * other.y) - (y * other.x);
  }

  [[nodiscard]] Real magnitude() const { return std::sqrt((x * x) + (y * y)); }

  [[nodiscard]] Vector2D normalized() const {
    Real mag = magnitude();
    return (mag > 0.0) ? Vector2D{x / mag, y / mag} : Vector2D{0.0, 0.0};
  }

  [[nodiscard]] Vector2D perpendicular() const {
    return {-y, x}; // 90-degree counter-clockwise rotation
  }
};

} // namespace bem::types
