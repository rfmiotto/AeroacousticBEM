#pragma once

#include <optional>
#include <unordered_map>

#include "domain/bem/DofMapper.hpp"
#include "domain/geometry/FlatPlate.hpp"
#include "domain/physics/BoundaryConditionProvider.hpp"
#include "foundation/types/BoundaryValueTypes.hpp"
#include "foundation/types/GeometryTypes.hpp"
#include "physics/RegionBCProviders.hpp"

namespace bem::domain::bem {

using domain::geometry::FlatPlate2D;
using foundation::dof::DofMapper;
using types::BCType;
using types::BoundaryCondition;
using types::ElementRegion;

/**
 * @brief Manages the assignment and enforcement of boundary conditions (BCs)
 * on the global DOFs of the BEM system.
 *
 * - BCs are assigned per element region (UPSTREAM, AIRFOIL, DOWNSTREAM).
 * - Nodes shared between different regions resolve conflicts by prioritizing
 *   AIRFOIL > DOWNSTREAM > UPSTREAM.
 * - Supports Dirichlet, Neumann, and Kutta boundary conditions.
 * - Uses BoundaryConditionProvider objects to handle region-specific BCs,
 *   allowing either fixed values or functions of coordinates (and possibly DOF
 *   context).
 */
class BoundaryConditionManager {
public:
  BoundaryConditionManager(const FlatPlate2D *geometry,
                           const DofMapper *mapper,
                           IBoundaryConditionProvider *upstream,
                           IBoundaryConditionProvider *airfoil,
                           IBoundaryConditionProvider *downstream);

  /// Assign BCs across all regions (idempotent)
  void assignAllRegionBCs();

  /// Update BCs only for the given region
  void updateRegionBC(ElementRegion region);

  /// Get BC for a global DOF (nullopt if not assigned)
  [[nodiscard]] std::optional<BoundaryCondition>
  getDofBC(std::size_t globalIdx) const;

  /// Get BC for an element's local DOF
  [[nodiscard]] std::optional<BoundaryCondition>
  getElementDofBC(std::size_t elemIdx, std::size_t localIdx) const;

  /// Get all global BCs (for solver assembly)
  [[nodiscard]] const std::unordered_map<std::size_t, BoundaryCondition> &
  allBCs() const {
    return dofBCs_;
  }

  /// Clear all assigned BCs
  void clear();

private:
  const FlatPlate2D *geometry_;
  const DofMapper *mapper_;

  // Stores final BCs per global DOF
  std::unordered_map<std::size_t, BoundaryCondition> dofBCs_;

  // Region providers
  IBoundaryConditionProvider *upstreamProvider_;
  IBoundaryConditionProvider *airfoilProvider_;
  IBoundaryConditionProvider *downstreamProvider_;

  void assignRegion(ElementRegion region);
  void resolveConflicts();
  static int regionPriority(ElementRegion region);
};

} // namespace bem::domain::bem
