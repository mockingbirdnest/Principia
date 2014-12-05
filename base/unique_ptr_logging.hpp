#pragma once

#include <memory>

// Logs the pointer.  No transfer of ownership.
template <typename T, typename U>
std::ostream& operator<<(std::ostream& out,
                         std::unique_ptr<T, U> const& pointer) {
  return out << pointer.get();
}
