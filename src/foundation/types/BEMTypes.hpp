#pragma once

#include "MathObjectTypes.hpp"

namespace bem::types {

// Element identification
using ElementId = std::size_t;
using NodeId = std::size_t;

// Boundary element types
enum class ElementType : std::uint8_t {
  LINEAR,    // 2-node linear element
  QUADRATIC, // 3-node quadratic element
  CUBIC      // 4-node cubic element
};

// Integration rule types
enum class IntegrationType : std::uint8_t {
  GAUSS,     // Gauss-Legendre quadrature
  ADAPTIVE,  // Adaptive quadrature
  ANALYTICAL // Analytical integration (when possible)
};

// Integration parameters
struct IntegrationParameters {
  IntegrationType type{IntegrationType::GAUSS};
  int order{8};              // Number of quadrature points
  Real tolerance{1e-12};     // For adaptive integration
  int max_subdivisions{100}; // For adaptive integration
  bool use_singularity_treatment{true};
  Real singularity_distance{1e-6}; // Distance threshold for singularity
};

// Quadrature point
struct QuadraturePoint {
  Real xi{0.0};     // Parameter coordinate [-1, 1]
  Real weight{0.0}; // Integration weight
  Point2D point;    // Physical coordinate

  QuadraturePoint() = default;
  QuadraturePoint(Real xi_val, Real weight_val, const Point2D &point_val)
      : xi(xi_val), weight(weight_val), point(point_val) {
  }
};

// Solver parameters
struct SolverParameters {
  Real tolerance{1e-12};
  int max_iterations{1000};
  bool use_preconditioning{true};
  std::string solver_type{"direct"}; // "direct", "gmres", "bicgstab"
};

// BEM system components
struct SystemMatrixEntry {
  ElementId row;
  ElementId col;
  Complex value;

  SystemMatrixEntry(ElementId r, ElementId c, Complex v)
      : row(r), col(c), value(v) {
  }
};

// Problem parameters
struct ProblemParameters {
  Real frequency{0.0};
  Real wave_number{0.0};
};

// Solution at a point
struct FieldPoint {
  Point2D location;
  Complex potential;
  Complex pressure;  // p = -iωρφ for acoustic problems
  Vector2D velocity; // v = ∇φ (as complex components)

  FieldPoint() = default;
  explicit FieldPoint(const Point2D &loc) : location(loc) {
  }
};

} // namespace bem::types
