#include "BoundaryConditionManager.hpp"
#include "foundation/exceptions/BEMException.hpp"

namespace bem::boundary {

using foundation::exceptions::BoundaryConditionException;

void BoundaryConditionManager::setDirichlet(ElementId id, const Complex &phi) {
  bc_map_[id] = {BoundaryConditionType::DIRICHLET,
                 BoundaryValue{.potential = phi}};
}

void BoundaryConditionManager::setNeumann(ElementId id,
                                          const Complex &dphi_dn) {
  bc_map_[id] = {BoundaryConditionType::NEUMANN,
                 BoundaryValue{.normal_derivative = dphi_dn}};
}

void BoundaryConditionManager::setRobin(ElementId id,
                                        const Complex &phi,
                                        const Complex &dphi_dn) {
  bc_map_[id] = {BoundaryConditionType::ROBIN,
                 BoundaryValue{.potential = phi, .normal_derivative = dphi_dn}};
}

void BoundaryConditionManager::clear(ElementId id) {
  bc_map_.erase(id);
}

bool BoundaryConditionManager::has(ElementId id) const {
  return bc_map_.contains(id);
}

BoundaryConditionType BoundaryConditionManager::getType(ElementId id) const {
  if (!has(id)) {
    throw BoundaryConditionException("No BC set for this element.");
  }
  return bc_map_.at(id).first;
}

const BoundaryValue &BoundaryConditionManager::getValue(ElementId id) const {
  if (!has(id)) {
    throw BoundaryConditionException("No BC set for this element.");
  }
  return bc_map_.at(id).second;
}

} // namespace bem::boundary
