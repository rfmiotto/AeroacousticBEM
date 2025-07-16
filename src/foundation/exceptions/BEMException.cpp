#include "BEMException.hpp"

namespace bem::foundation::exceptions {

BEMException::BEMException(const std::string &message)
    : std::runtime_error(message) {
}

BEMIntegrationException::BEMIntegrationException(const std::string &msg)
    : BEMException("Integration error: " + msg) {
}

BEMSingularityException::BEMSingularityException(const std::string &msg)
    : BEMException("Singularity encountered: " + msg) {
}

BEMConvergenceException::BEMConvergenceException(const std::string &msg)
    : BEMException("Convergence failure: " + msg) {
}

BoundaryConditionException::BoundaryConditionException(const std::string &msg)
    : BEMException("Boundary condition error: " + msg) {
}

} // namespace bem::foundation::exceptions
