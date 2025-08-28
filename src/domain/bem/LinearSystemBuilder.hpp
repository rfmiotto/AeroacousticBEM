#pragma once

#include "domain/physics/BoundaryConditionManager.hpp"

namespace bem::domain::assembly {

struct LinearSystem {
  Eigen::MatrixXcd A;              // Final coefficient matrix
  Eigen::VectorXcd b;              // RHS vector
  std::vector<int> unknownIndices; // DOFs being solved for
  std::vector<int> knownIndices;   // DOFs prescribed
};

/**
 * @brief Build the reduced linear system A x = b
 *        from assembled H, D and boundary conditions.
 *
 *        The discrete boundary integral equation is given by
 *        [H]{u} = [D]{q}
 *        H comes from the double-layer potential, while D comes
 *        from the single-layer potential.
 *        On a given boundary, we don’t know both u and q everywhere.
 *        The boundary condition manager decides which is known and
 *        which is unknown at each element (it splits the DOFs into
 *        knowns and unknowns)
 *        Then we rearange the system to solve for the unknown side.
 *
 * Convention:
 *   - Dirichlet BCs: known potential u, unknown q.
 *   - Neumann BCs: known q, unknown u.
 */
inline LinearSystem
buildLinearSystem(const Eigen::MatrixXcd &H,
                  const Eigen::MatrixXcd &D,
                  const bem::boundary::BoundaryConditionManager &bcMgr,
                  int nGlobalDofs) {
  Eigen::VectorXcd u = Eigen::VectorXcd::Zero(nGlobalDofs);
  Eigen::VectorXcd q = Eigen::VectorXcd::Zero(nGlobalDofs);

  std::vector<int> unknowns;
  std::vector<int> knowns;

  // Fill u, q with BC values where available
  for (int i = 0; i < nGlobalDofs; ++i) {
    if (bcMgr.has(i)) {
      auto type = bcMgr.getType(i);
      const auto &val = bcMgr.getValue(i);

      if (type == bem::types::BoundaryConditionType::DIRICHLET) {
        u(i) = val.potential;
        unknowns.push_back(i); // q is unknown
      } else if (type == bem::types::BoundaryConditionType::NEUMANN) {
        q(i) = val.normal_derivative;
        unknowns.push_back(i); // u is unknown
      } else if (type == bem::types::BoundaryConditionType::ROBIN) {
        throw std::runtime_error("Robin BC not implemented in builder yet.");
      }
    } else {
      throw std::runtime_error("No BC set for DOF " + std::to_string(i));
    }
  }

  // Form residual: H u - D q = 0
  Eigen::VectorXcd rhs = D * q - H * u;

  // Build reduced matrix (subselecting unknown columns)
  Eigen::MatrixXcd A(H.rows(), unknowns.size());
  for (int j = 0; j < (int)unknowns.size(); ++j) {
    int idx = unknowns[j];
    A.col(j) =
        H.col(idx) *
            (bcMgr.getType(idx) == bem::types::BoundaryConditionType::NEUMANN
                 ? 1.0
                 : 0.0) -
        D.col(idx) *
            (bcMgr.getType(idx) == bem::types::BoundaryConditionType::DIRICHLET
                 ? 1.0
                 : 0.0);
  }

  LinearSystem sys;
  sys.A = A;
  sys.b = rhs;
  sys.unknownIndices = unknowns;
  sys.knownIndices = knowns;
  return sys;
}

} // namespace bem::domain::assembly
