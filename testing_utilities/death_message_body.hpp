#pragma once

#include <string>

#include "testing_utilities/death_message.hpp"

namespace principia {
namespace testing_utilities {

inline std::string DeathMessage(std::string const& s) {
#ifdef _DEBUG
  return "";
#else
  return s;
#endif
}

}  // namespace testing_utilities
}  // namespace principia
