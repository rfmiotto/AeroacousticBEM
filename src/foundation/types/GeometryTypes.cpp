#include "foundation/types/GeometryTypes.hpp"
#include "foundation/utils/MathUtils.hpp"

namespace bem::types {
using namespace bem::foundation::utils::MathUtils;

void Element::setElementRegion(ElementRegion r) {
  region = r;
}

ElementRegion Element::getElementRegion() const {
  return region;
}

Real Element::length() const {
  return distance(start, end);
}

Point2D Element::midpoint() const {
  return Point2D{0.5 * (start.x + end.x), 0.5 * (start.y + end.y)};
}

Vector2D Element::tangent() const {
  return Vector2D{end.x - start.x, end.y - start.y};
}

Vector2D Element::normal() const {
  Vector2D t = tangent();
  return Vector2D{-t.y, t.x};
}

// Parametric point with xi ∈ [-1, 1]
Point2D Element::parametricPoint(Real xi) const {
  const Real a = 0.5 * (1.0 - xi);
  const Real b = 0.5 * (1.0 + xi);
  return Point2D{(a * start.x) + (b * end.x), (a * start.y) + (b * end.y)};
}

// Jacobian for numerical integration in the [-1, 1] parametric domain
Real Element::jacobian() const {
  return 0.5 * length();
}

bool Element::contains(const Point2D &pt, Real tol) const {
  // Element vector
  const Real dx = end.x - start.x;
  const Real dy = end.y - start.y;

  // vector from element start to point
  const Real dxp = pt.x - start.x;
  const Real dyp = pt.y - start.y;

  // If the cross product between these two vectors is not zero,
  // it is not aligned → point is not on the line segment
  const Real cross = (dx * dyp) - (dy * dxp);

  if (std::abs(cross) > tol) {
    return false;
  }

  // Dot produt to check if point is within element extremes
  const Real dot = (dx * dxp) + (dy * dyp);
  const Real len2 = (dx * dx) + (dy * dy);

  return (dot >= -tol) && (dot <= len2 + tol);
}

std::vector<Real> Element::shapeFunction(Real xi) const {
  switch (type) {
  case ElementType::CONSTANT:
    return {1.0};

  case ElementType::LINEAR: {
    // Nodes at xi = -1, +1
    Real N1 = 0.5 * (1.0 - xi);
    Real N2 = 0.5 * (1.0 + xi);
    return {
        N1,
        N2,
    };
  }

  case ElementType::QUADRATIC: {
    // Nodes at xi = -1, 0, +1
    Real N1 = 0.5 * xi * (xi - 1.0);
    Real N2 = 1.0 - (xi * xi);
    Real N3 = 0.5 * xi * (xi + 1.0);
    return {
        N1,
        N2,
        N3,
    };
  }

  case ElementType::CUBIC: {
    // Lagrange interpolation through xi = -1, -1/3, 1/3, 1
    constexpr Real X1 = -1.0;
    constexpr Real X2 = -1.0 / 3.0;
    constexpr Real X3 = 1.0 / 3.0;
    constexpr Real X4 = 1.0;
    auto L = [&](Real xi, Real x_i) {
      Real xs[4] = {X1, X2, X3, X4};
      Real val = 1.0;
      for (double x : xs) {
        if (x != x_i) {
          val *= (xi - x) / (x_i - x);
        }
      }
      return val;
    };
    return {
        L(xi, X1),
        L(xi, X2),
        L(xi, X3),
        L(xi, X4),
    };
  }
  }
  return {};
}

std::vector<Real> Element::shapeFunctionDerivative(Real xi) const {
  switch (type) {
  case ElementType::CONSTANT:
    return {0.0};

  case ElementType::LINEAR:
    return {-0.5, 0.5};

  case ElementType::QUADRATIC:
    return {xi - 0.5, -2.0 * xi, xi + 0.5};

  case ElementType::CUBIC: {
    constexpr Real X1 = -1.0;
    constexpr Real X2 = -1.0 / 3.0;
    constexpr Real X3 = 1.0 / 3.0;
    constexpr Real X4 = 1.0;
    auto dL = [&](Real xi, int i) {
      Real xs[4] = {X1, X2, X3, X4};
      Real sum = 0.0;
      for (int m = 0; m < 4; ++m) {
        if (m == i) {
          continue;
        }
        Real term = 1.0 / (xs[i] - xs[m]);
        for (int j = 0; j < 4; ++j) {
          if (j == i || j == m) {
            continue;
          }
          term *= (xi - xs[j]) / (xs[i] - xs[j]);
        }
        sum += term;
      }
      return sum;
    };
    return {dL(xi, 0), dL(xi, 1), dL(xi, 2), dL(xi, 3)};
  }
  }
  return {};
}

} // namespace bem::types
