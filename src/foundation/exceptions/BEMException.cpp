#include "BEMException.hpp"

namespace bem::foundation::exceptions {

BEMException::BEMException(const std::string &message)
    : std::runtime_error(message) {}

} // namespace bem::foundation::exceptions
