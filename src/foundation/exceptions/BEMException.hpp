#pragma once

#include <stdexcept>
#include <string>

namespace bem::foundation::exceptions {

class BEMException : public std::runtime_error {
public:
  explicit BEMException(const std::string &message);
};

class BEMIntegrationException : public BEMException {
public:
  explicit BEMIntegrationException(const std::string &msg);
};

class BEMSingularityException : public BEMException {
public:
  explicit BEMSingularityException(const std::string &msg);
};

class BEMConvergenceException : public BEMException {
public:
  explicit BEMConvergenceException(const std::string &msg);
};

} // namespace bem::foundation::exceptions
