#pragma once

#include "base/monostable.hpp"

namespace principia {
namespace base {
namespace _monostable {
namespace internal {

inline void Monostable::Flop() {
  transient_ = false;
}

inline Monostable::operator bool() const {
  return transient_;
}

}  // namespace internal
}  // namespace _monostable
}  // namespace base
}  // namespace principia
