#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "foundation/types/GeometryTypes.hpp"

namespace bem::foundation::dof {

using bem::types::Element;
using bem::types::ElementType;
using bem::types::Point2D;
using bem::types::Real;

/**
 * @brief DofMapper maps element-local DOF indices to global DOF indices and
 * vice-versa for a structured/sequential list of boundary elements. This
 * implementation takes advantage of the fact that all elements share the same
 * ElementType and are stored in sequential order around the geometry.
 *
 * Key design choices:
 * - Linear pass over elements: uses the sequential ordering to avoid costly
 * global-coordinate hashing. For element i>0 we assume the first local
 * node is shared with the previous element's last local node (for
 * LINEAR/QUADRATIC/CUBIC). For CONSTANT elements no node sharing occurs.
 * - Precomputed element->global index table: vector<vector<std::size_t>>
 * elementDofs_ so assembly of global matrices is cheap and cache friendly.
 * - Global->(element,local) reverse mapping for quick lookups of ownership.
 * - Simple API for assembling matrices: getElementGlobalIds(elementIdx)
 * returns a span (vector reference) of global DOF indices for that element.
 */
class DofMapper {
public:
  explicit DofMapper(const std::vector<Element> &elements);

  // Build the mapping. Separated from constructor for clarity and possible
  // re-builds when BCs change.
  void build();

  // Get number of global DOFs
  [[nodiscard]] std::size_t nGlobalDofs() const {
    return globalPoints_.size();
  }

  // Get a const reference to the global unique DOF coordinates
  [[nodiscard]] const std::vector<Point2D> &globalPoints() const {
    return globalPoints_;
  }

  // For a given element index, returns the vector of global DOF indices in
  // the same local ordering produced by Element::getNodes().
  [[nodiscard]] const std::vector<std::size_t> &
  getElementGlobalIds(std::size_t elemIdx) const;

  // For a global DOF, returns the list of (elementIdx, localIdx) that reference
  // it.
  [[nodiscard]] const std::vector<std::pair<std::size_t, std::size_t>> &
  getGlobalOwners(std::size_t globalIdx) const;

  // Convenience: return global index for pair (elementIdx, localIdx)
  [[nodiscard]] std::size_t elementLocalToGlobal(std::size_t elementIdx,
                                                 std::size_t localIdx) const;

  // Return list of global DOF indices that have Dirichlet/essential BCs.
  // The caller is responsible for filling this vector (e.g. from boundary
  // region information or external BC specification). We store it here only
  // for convenience during assembly.
  [[nodiscard]] const std::vector<std::size_t> &dirichletIndices() const {
    return dirichletIndices_;
  }

private:
  const std::vector<Element> *elementsPtr_;
  ElementType elementType_; // will be inferred from first element

  // Unique global DOF coordinates
  std::vector<Point2D> globalPoints_;

  // maps from each element → its list of global DOF indices (local node order
  // preserved) “For element e, what are the global indices of its local nodes?”
  std::vector<std::vector<std::size_t>> elementDofs_;

  // globalId -> list of (elementId, localNodeIdx)
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> globalOwners_;

  // Indices with Dirichlet BCs
  std::vector<std::size_t> dirichletIndices_;

  // Add all nodes from the first element to initialize the mapping
  void addFirstElementNodes();

  // Add nodes from subsequent elements, handling shared nodes properly
  void addSubsequentElementNodes(std::size_t elemIdx);

  // Ensure closure of the domain by merging first and last node if they are
  // equal
  void handleClosure();

  // Build reverse mapping from global DOF to element-local ownership
  void buildReverseMapping();

  // Reset internal storage before rebuilding
  void reset();
};

} // namespace bem::foundation::dof
