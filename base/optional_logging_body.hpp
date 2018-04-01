
#pragma once

#include "base/optional_logging.hpp"

template<typename T>
std::ostream& operator<<(std::ostream& out,
                         std::optional<T> const& optional) {
  if (optional) {
    return out << *optional;
  } else {
    return out << "nullopt";
  }
}
