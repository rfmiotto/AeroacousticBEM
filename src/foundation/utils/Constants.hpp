#pragma once

#include "foundation/types/MathObjectTypes.hpp"
#include <numbers>

namespace bem::foundation::utils::Constants {

using types::Complex;
using types::Real;

// Mathematical constants
static constexpr Real PI = std::numbers::pi;
static constexpr Real E = std::numbers::e;
static constexpr Real EULER_MASCHERONI = std::numbers::egamma;
static constexpr Real SQRT_PI = 1.7724538509055159;
static constexpr Real INV_PI = std::numbers::inv_pi;
static constexpr Real INV_2PI = std::numbers::inv_pi / 2.0;
static constexpr Real INV_4PI = std::numbers::inv_pi / 4.0;

// Complex constants
static constexpr Complex I{0.0, 1.0};             // Imaginary unit
static constexpr Complex I_4{0.0, 0.25};          // i/4
static constexpr Complex INV_4PI_I{0.0, INV_4PI}; // i/(4π)

// Numerical constants
static constexpr Real MACHINE_EPSILON = std::numeric_limits<Real>::epsilon();
static constexpr Real SQRT_MACHINE_EPSILON = 1.4901161193847656e-08;
static constexpr Real INTEGRATION_TOLERANCE = 1e-12;
static constexpr Real SOLVER_TOLERANCE = 1e-12;
static constexpr Real GEOMETRY_TOLERANCE = 1e-14;

// BEM-specific constants
static constexpr Real SINGULARITY_DISTANCE_THRESHOLD = 1e-6;
static constexpr Real FAR_FIELD_DISTANCE_FACTOR = 10.0;
static constexpr int DEFAULT_INTEGRATION_ORDER = 4;
static constexpr int DEFAULT_MAX_SUBDIVISIONS = 100;
static constexpr bool USE_SINGULARITY_TREATMENT = true;

// Hankel function computation constants
static constexpr Real HANKEL_LARGE_ARG_THRESHOLD = 8.0;
static constexpr Real HANKEL_SMALL_ARG_THRESHOLD = 1e-8;
static constexpr int HANKEL_SERIES_MAX_TERMS = 100;

} // namespace bem::foundation::utils::Constants
