#pragma once

#include "base/unique_ptr_logging.hpp"

#include <memory>

template<typename T, typename U>
std::ostream& operator<<(std::ostream& out,
                         std::unique_ptr<T, U> const& pointer) {
  return out << pointer.get();
}
