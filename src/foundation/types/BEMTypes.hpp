#pragma once

#include "ComplexTypes.hpp"
#include "GeometryTypes.hpp"
#include <functional>

namespace bem::types {

// BEM system components
struct SystemMatrixEntry {
  ElementId row;
  ElementId col;
  Complex value;

  SystemMatrixEntry(ElementId r, ElementId c, Complex v)
      : row(r), col(c), value(v) {}
};

// Boundary value
struct BoundaryValue {
  Complex potential;         // φ (acoustic potential)
  Complex normal_derivative; // ∂φ/∂n (normal velocity)
  bool potential_known{false};
  bool derivative_known{false};

  BoundaryValue() = default;
  explicit BoundaryValue(Complex pot, bool pot_known = true)
      : potential(pot), potential_known(pot_known) {}
};

// Solution at a point
struct FieldPoint {
  Point2D location;
  Complex potential;
  Complex pressure;  // p = -iωρφ for acoustic problems
  Vector2D velocity; // v = ∇φ (as complex components)

  FieldPoint() = default;
  explicit FieldPoint(const Point2D &loc) : location(loc) {}
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

// Solver parameters
struct SolverParameters {
  Real tolerance{1e-12};
  int max_iterations{1000};
  bool use_preconditioning{true};
  std::string solver_type{"direct"}; // "direct", "gmres", "bicgstab"
};

// Problem parameters
struct ProblemParameters {
  Frequency frequency{0.0};
  WaveNumber wave_number{0.0};
  Real density{1.225};     // Air density (kg/m³)
  Real sound_speed{343.0}; // Sound speed (m/s)
};

// Function types for boundary conditions
using BoundaryFunction = std::function<Complex(const Point2D &)>;
using NormalFunction = std::function<Vector2D(const Point2D &)>;

} // namespace bem::types
