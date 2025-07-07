#pragma once

#include "foundation/types/MathObjectTypes.hpp"
#include <numbers>

namespace bem::foundation::utils::Constants {

// Mathematical constants
static constexpr types::Real PI = std::numbers::pi;
static constexpr types::Real E = std::numbers::e;
static constexpr types::Real EULER_MASCHERONI = std::numbers::egamma;
static constexpr types::Real SQRT_PI = 1.7724538509055159;
static constexpr types::Real INV_PI = std::numbers::inv_pi;
static constexpr types::Real INV_2PI = std::numbers::inv_pi / 2.0;
static constexpr types::Real INV_4PI = std::numbers::inv_pi / 4.0;

// Complex constants
static constexpr types::Complex I{0.0, 1.0};             // Imaginary unit
static constexpr types::Complex I_4{0.0, 0.25};          // i/4
static constexpr types::Complex INV_4PI_I{0.0, INV_4PI}; // i/(4π)

// Numerical constants
static constexpr types::Real MACHINE_EPSILON =
    std::numeric_limits<types::Real>::epsilon();
static constexpr types::Real SQRT_MACHINE_EPSILON = 1.4901161193847656e-08;
static constexpr types::Real INTEGRATION_TOLERANCE = 1e-12;
static constexpr types::Real SOLVER_TOLERANCE = 1e-12;
static constexpr types::Real GEOMETRY_TOLERANCE = 1e-14;

// BEM-specific constants
static constexpr types::Real SINGULARITY_DISTANCE_THRESHOLD = 1e-6;
static constexpr types::Real FAR_FIELD_DISTANCE_FACTOR = 10.0;
static constexpr int DEFAULT_INTEGRATION_ORDER = 8;
static constexpr int MAX_INTEGRATION_ORDER = 32;

// Hankel function computation constants
static constexpr types::Real HANKEL_LARGE_ARG_THRESHOLD = 8.0;
static constexpr types::Real HANKEL_SMALL_ARG_THRESHOLD = 1e-8;
static constexpr int HANKEL_SERIES_MAX_TERMS = 100;

} // namespace bem::foundation::utils::Constants
