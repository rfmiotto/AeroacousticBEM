#include "BoundaryConditionManager.hpp"
#include <stdexcept>

namespace bem::domain::bem {

BoundaryConditionManager::BoundaryConditionManager(
    const FlatPlate2D *geometry,
    const DofMapper *mapper,
    IBoundaryConditionProvider *upstream,
    IBoundaryConditionProvider *airfoil,
    IBoundaryConditionProvider *downstream)
    : geometry_(geometry), mapper_(mapper), upstreamProvider_(upstream),
      airfoilProvider_(airfoil), downstreamProvider_(downstream) {
}

void BoundaryConditionManager::assignAllRegionBCs() {
  dofBCs_.clear();
  assignRegion(ElementRegion::UPSTREAM);
  assignRegion(ElementRegion::AIRFOIL);
  assignRegion(ElementRegion::DOWNSTREAM);
  resolveConflicts();
}

void BoundaryConditionManager::updateRegionBC(ElementRegion region) {
  assignRegion(region);
  resolveConflicts();
}

void BoundaryConditionManager::assignRegion(ElementRegion region) {
  const auto &elements = geometry_->getElements();

  IBoundaryConditionProvider *provider = nullptr;
  switch (region) {
  case ElementRegion::UPSTREAM:
    provider = upstreamProvider_;
    break;
  case ElementRegion::AIRFOIL:
    provider = airfoilProvider_;
    break;
  case ElementRegion::DOWNSTREAM:
    provider = downstreamProvider_;
    break;
  }

  if (provider == nullptr) {
    return;
  }

  for (std::size_t e = 0; e < elements.size(); ++e) {
    if (elements[e].getElementRegion() != region) {
      continue;
    }

    for (std::size_t g : mapper_->getElementGlobalIds(e)) {
      const auto &pt = mapper_->globalPoints()[g];
      dofBCs_[g] = provider->evaluate(pt.x); // overwrite or insert
    }
  }
}

std::optional<BoundaryCondition>
BoundaryConditionManager::getDofBC(std::size_t globalIdx) const {
  auto it = dofBCs_.find(globalIdx);
  if (it != dofBCs_.end()) {
    return it->second;
  }
  return std::nullopt;
}

std::optional<BoundaryCondition>
BoundaryConditionManager::getElementDofBC(std::size_t elemIdx,
                                          std::size_t localIdx) const {
  const auto &globalIds = mapper_->getElementGlobalIds(elemIdx);
  if (localIdx >= globalIds.size()) {
    throw std::out_of_range("Invalid local DOF index in getElementDofBC");
  }
  std::size_t g = globalIds[localIdx];
  return getDofBC(g);
}

void BoundaryConditionManager::clear() {
  dofBCs_.clear();
}

void BoundaryConditionManager::resolveConflicts() {
  const auto &elements = geometry_->getElements();

  for (std::size_t g = 0; g < mapper_->nGlobalDofs(); ++g) {
    const auto &owners = mapper_->getGlobalOwners(g);
    if (owners.size() < 2) {
      continue;
    }

    // Select region with max priority
    ElementRegion chosenRegion = ElementRegion::UPSTREAM;
    int maxPriority = -1;
    for (const auto &[elemIdx, localIdx] : owners) {
      ElementRegion r = elements[elemIdx].getElementRegion();
      int p = regionPriority(r);
      if (p > maxPriority) {
        maxPriority = p;
        chosenRegion = r;
      }
    }

    IBoundaryConditionProvider *provider = nullptr;
    switch (chosenRegion) {
    case ElementRegion::UPSTREAM:
      provider = upstreamProvider_;
      break;
    case ElementRegion::AIRFOIL:
      provider = airfoilProvider_;
      break;
    case ElementRegion::DOWNSTREAM:
      provider = downstreamProvider_;
      break;
    }

    if (provider != nullptr) {
      const auto &pt = mapper_->globalPoints()[g];
      dofBCs_[g] = provider->evaluate(pt.x);
    }
  }
}

int BoundaryConditionManager::regionPriority(ElementRegion region) {
  switch (region) {
  case ElementRegion::AIRFOIL:
    return 2;
  case ElementRegion::DOWNSTREAM:
    return 1;
  case ElementRegion::UPSTREAM:
  default:
    return 0;
  }
}

} // namespace bem::domain::bem
