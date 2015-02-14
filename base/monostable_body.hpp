#pragma once

#include "base/monostable.hpp"

namespace principia {
namespace base {

inline void Monostable::Flop() {
  transient_ = false;
}

inline Monostable::operator bool() const {
  return transient_;
}

}  // namespace base
}  // namespace principia
