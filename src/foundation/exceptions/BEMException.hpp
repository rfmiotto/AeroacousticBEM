#pragma once

#include <stdexcept>
#include <string>

namespace bem::foundation::exceptions {

class BEMException : public std::runtime_error {
public:
  explicit BEMException(const std::string &message);
};

} // namespace bem::foundation::exceptions
