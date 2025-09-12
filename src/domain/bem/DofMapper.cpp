#include "DofMapper.hpp"
#include "foundation/utils/MathUtils.hpp"

namespace bem::foundation::dof {

using bem::foundation::utils::MathUtils::pointsEqual;

DofMapper::DofMapper(const std::vector<Element> &elements)
    : elementsPtr_(&elements) {

  if (!elements.empty()) {
    elementType_ = elements.front().type;
  }

  build();
}

void DofMapper::addFirstElementNodes() {
  const auto nodes = (*elementsPtr_)[0].getNodes();
  for (const auto &node : nodes) {
    globalPoints_.push_back(node);
    elementDofs_[0].push_back(globalPoints_.size() - 1);
  }
}

// From the second element on, use sequential-sharing assumption:
// - For CONSTANT elements: all nodes are unique (no sharing)
// - For other element types: assume first local node of element i equals
//   last local node of element i-1 (so we skip adding the first node and
//   reuse the previous last global index)
void DofMapper::addSubsequentElementNodes(std::size_t e) {
  const auto nodes = (*elementsPtr_)[e].getNodes();

  // For constant elements, add all nodes directly: each midpoint is a unique
  // DOF
  if (elementType_ == ElementType::CONSTANT) {
    for (const auto &node : nodes) {
      globalPoints_.push_back(node);
      elementDofs_[e].push_back(globalPoints_.size() - 1);
    }
    return;
  }

  // For higher-order elements, first node matches previous element's last node
  const std::size_t prevLastGlobal = elementDofs_[e - 1].back();
  elementDofs_[e].push_back(prevLastGlobal);

  // Add the remaining nodes as new global points
  for (std::size_t local = 1; local < nodes.size(); ++local) {
    globalPoints_.push_back(nodes[local]);
    elementDofs_[e].push_back(globalPoints_.size() - 1);
  }
}

// Circular closure: if last element's last node duplicates first element's
// first node, merge them. This is common when domain is closed.
void DofMapper::handleClosure() {
  if (globalPoints_.size() < 2) {
    return;
  }

  const Point2D &firstPt = globalPoints_.front();
  const Point2D &lastPt = globalPoints_.back();

  // If first and last points are equal, merge them into one global DOF
  if (pointsEqual(firstPt, lastPt)) {
    const std::size_t lastIdx = globalPoints_.size() - 1;
    for (auto &elemVec : elementDofs_) {
      for (auto &g : elemVec) {
        if (g == lastIdx) {
          g = 0;
        }
      }
    }
    globalPoints_.pop_back();
  }
}

// Build reverse mapping: global -> list of owners
void DofMapper::buildReverseMapping() {
  globalOwners_.resize(globalPoints_.size());
  for (std::size_t e = 0; e < elementsPtr_->size(); ++e) {
    for (std::size_t local = 0; local < elementDofs_[e].size(); ++local) {
      const std::size_t g = elementDofs_[e][local];
      if (g >= globalOwners_.size()) {
        throw std::runtime_error("Internal error building globalOwners_");
      }
      globalOwners_[g].emplace_back(e, local);
    }
  }
}

void DofMapper::reset() {
  elementDofs_.clear();
  elementDofs_.resize(elementsPtr_->size());
  globalPoints_.clear();
  globalOwners_.clear();
}

void DofMapper::build() {
  reset();
  if (elementsPtr_->empty()) {
    return;
  }

  addFirstElementNodes();

  for (std::size_t e = 1; e < elementsPtr_->size(); ++e) {
    addSubsequentElementNodes(e);
  }

  handleClosure();

  buildReverseMapping();
}

const std::vector<std::size_t> &
DofMapper::getElementGlobalIds(std::size_t elemIdx) const {
  if (elemIdx >= elementDofs_.size()) {
    throw std::out_of_range(
        "Element index out of range in getElementGlobalIds");
  }
  return elementDofs_[elemIdx];
}

const std::vector<std::pair<std::size_t, std::size_t>> &
DofMapper::getGlobalOwners(std::size_t globalIdx) const {
  if (globalIdx >= globalOwners_.size()) {
    throw std::out_of_range("Global index out of range in getGlobalOwners");
  }
  return globalOwners_[globalIdx];
}

std::size_t DofMapper::elementLocalToGlobal(std::size_t elementIdx,
                                            std::size_t localIdx) const {
  if (elementIdx >= elementDofs_.size()) {
    throw std::out_of_range(
        "Element index out of range in elementLocalToGlobal");
  }
  const auto &vec = elementDofs_[elementIdx];
  if (localIdx >= vec.size()) {
    throw std::out_of_range("Local index out of range in elementLocalToGlobal");
  }
  return vec[localIdx];
}

} // namespace bem::foundation::dof
