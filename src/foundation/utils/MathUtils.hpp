#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <cmath>
#include <complex>
#include <utility> // For std::pair
#include <vector>

#include "foundation/exceptions/BEMException.hpp"
#include "foundation/types/ComplexTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/Constants.hpp"

namespace bem::foundation::utils::MathUtils {

// Function declarations
types::Real distance(const types::Point2D &p1, const types::Point2D &p2);
types::Real distanceSquared(const types::Point2D &p1, const types::Point2D &p2);
types::Real angle(const types::Point2D &p1, const types::Point2D &p2);
bool isZero(double value, double tol = Constants::MACHINE_EPSILON);
bool isEqual(double v1, double v2, double tol = Constants::MACHINE_EPSILON);
types::Real magnitude(types::Complex z);
types::Real phase(types::Complex z);
types::Complex hankelFirstKind(int order, double z);
types::Complex greensFunction2D(const types::Point2D &source,
                                const types::Point2D &field, double k);
types::Complex greensFunctionDerivative2D(const types::Point2D &source,
                                          const types::Point2D &field,
                                          const types::Point2D &normal,
                                          double k);
std::pair<std::vector<double>, std::vector<double>>
gaussLegendreQuadrature(int n_points);
std::pair<types::Point2D, double>
lineIntegrationTransform(const types::Point2D &start, const types::Point2D &end,
                         double xi);

} // namespace bem::foundation::utils::MathUtils

#endif // MATHUTILS_HPP
