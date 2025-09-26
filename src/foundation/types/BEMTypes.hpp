#pragma once

#include "MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"

namespace bem::types {

using types::Complex;
using types::Real;
namespace Constants = bem::foundation::utils::Constants;

// Element identification
using ElementId = std::size_t;
using NodeId = std::size_t;

// Integration rule types
enum class IntegrationType : std::uint8_t {
  GAUSS,    // Gauss-Legendre quadrature
  ADAPTIVE, // Adaptive quadrature
};

// Integration parameters
struct IntegrationParameters {
  IntegrationType type{IntegrationType::GAUSS};
  int order{
      Constants::DEFAULT_INTEGRATION_ORDER}; // Number of quadrature points
  Real tolerance{Constants::INTEGRATION_TOLERANCE};          // For adaptive
                                                             // integration
  int max_subdivisions{Constants::DEFAULT_MAX_SUBDIVISIONS}; // For adaptive
                                                             // integration
  bool use_singularity_treatment{Constants::USE_SINGULARITY_TREATMENT};
  Real singularity_distance{
      Constants::SINGULARITY_DISTANCE_THRESHOLD}; // Distance threshold for
                                                  // singularity
};

// Solver parameters
struct SolverParameters {
  Real tolerance{1e-12};
  int max_iterations{1000};
  bool use_preconditioning{true};
  std::string solver_type{"direct"}; // "direct", "gmres", "bicgstab"
};

} // namespace bem::types
