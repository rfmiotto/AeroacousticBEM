#pragma once

#include <cmath>

#include "foundation/types/MathObjectTypes.hpp"
#include "foundation/utils/Constants.hpp"

using namespace bem::foundation::utils::Constants;

namespace bem::foundation::utils::MathUtils {

// Function declarations
types::Real distance(const types::Point2D &p1, const types::Point2D &p2);
types::Real distanceSquared(const types::Point2D &p1, const types::Point2D &p2);
types::Real angle(const types::Point2D &p1, const types::Point2D &p2);
bool isZero(double value, double tol = MACHINE_EPSILON);
bool isEqual(double v1, double v2, double tol = MACHINE_EPSILON);
types::Real magnitude(types::Complex z);
types::Real phase(types::Complex z);

} // namespace bem::foundation::utils::MathUtils
