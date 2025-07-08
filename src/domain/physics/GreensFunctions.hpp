#pragma once

#include "foundation/types/MathObjectTypes.hpp"

namespace bem::domain::physics {

using bem::types::Complex;
using bem::types::Point2D;
using bem::types::Real;
using bem::types::Vector2D;

/**
 * @brief Computes the Hankel function of the first kind.
 *
 * This is typically used for representing outgoing wave solutions in 2D
 * Helmholtz problems.
 *
 * @param order Order of the Hankel function (usually 0 or 1)
 * @param z     Complex argument
 * @ return Value of the Hankel function of the first kind
 */
Complex hankelFirstKind(int order, Real z);

/**
 * @brief Computes the 2D free-space Green’s function for the Helmholtz
 * equation.
 * This represents the fundamental solution of the Helmholtz equation in
 * two dimensions, satisfying outgoing wave conditions.
 *
 * Free-field Green's function G(x,y) for 2D Helmholtz equation
 * G(x,y) = (i/4) * H_0^(1)(k|x-y|)
 * where H_0^(1) is the Hankel function of first kind, order 0
 *
 * @param source    Location of the source point
 * @param observer  Location of the observer (evaluation) point
 * @param k         Wave number (must be positive)
 * @return          Value of the Green’s function at the observer point
 */
Complex
greensFunction2D(const Point2D &source, const Point2D &observer, Real k);

/**
 * @brief Computes the normal derivative of the 2D Green’s function.
 *
 * This is used for boundary integral formulations, such as the Neumann
 * boundary condition in Helmholtz problems.
 *
 * Normal derivative of Green's function
 * ∂G/∂n = (i/4) * k * H_1^(1)(k|x-y|) * (∂r/∂n)
 *
 * @param source    Location of the source point
 * @param observer  Location of the observer point (evaluation point)
 * @param normal    Outward unit normal vector at the source point
 * @param k         Wave number
 * @return          Normal derivative of the Green’s function
 */
Complex greensFunctionNormalDerivative2D(const Point2D &source,
                                         const Point2D &observer,
                                         const Vector2D &normal,
                                         Real k);

} // namespace bem::domain::physics
