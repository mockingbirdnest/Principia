#include "base/monostable.hpp"

namespace principia {
namespace base {

void Monostable::Flop() {
  transient_ = false;
}

Monostable::operator bool() const {
  return transient_;
}

}  // namespace base
}  // namespace principia
