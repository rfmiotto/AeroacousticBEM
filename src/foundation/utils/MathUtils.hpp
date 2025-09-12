#pragma once

#include <cmath>

#include "foundation/types/MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"

namespace bem::foundation::utils::MathUtils {

using bem::foundation::utils::Constants::GEOMETRY_TOLERANCE;
using bem::foundation::utils::Constants::MACHINE_EPSILON;
using types::Complex;
using types::Point2D;
using types::Real;

// Function declarations
Real distance(const Point2D &p1, const Point2D &p2);
Real distanceSquared(const Point2D &p1, const Point2D &p2);
Real angle(const Point2D &p1, const Point2D &p2);
bool isZero(double value, double tol = MACHINE_EPSILON);
bool isEqual(double v1, double v2, double tol = MACHINE_EPSILON);
bool pointsEqual(const Point2D &p1,
                 const Point2D &p2,
                 double tol = GEOMETRY_TOLERANCE);
Real magnitude(Complex z);
Real phase(Complex z);

} // namespace bem::foundation::utils::MathUtils
