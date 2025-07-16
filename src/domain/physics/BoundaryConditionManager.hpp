#pragma once

#include "foundation/types/BEMTypes.hpp"
#include "foundation/types/BoundaryValueTypes.hpp"
#include "foundation/types/MathObjectTypes.hpp"

namespace bem::boundary {

using bem::types::BoundaryConditionType;
using bem::types::BoundaryValue;
using bem::types::Complex;
using bem::types::ElementId;

/**
 * @brief Manages boundary conditions for a set of elements.
 */
class BoundaryConditionManager {
public:
  void setDirichlet(ElementId id, const Complex &phi);
  void setNeumann(ElementId id, const Complex &dphi_dn);
  void setRobin(ElementId id, const Complex &phi, const Complex &dphi_dn);
  void clear(ElementId id);

  [[nodiscard]] bool has(ElementId id) const;

  [[nodiscard]] BoundaryConditionType getType(ElementId id) const;
  [[nodiscard]] const BoundaryValue &getValue(ElementId id) const;

private:
  std::unordered_map<ElementId, std::pair<BoundaryConditionType, BoundaryValue>>
      bc_map_;
};

} // namespace bem::boundary
