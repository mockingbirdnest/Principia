#include "ksp_plugin/monostable.hpp"

namespace principia {
namespace ksp_plugin {

void Monostable::Flop() {
  transient_ = false;
}

Monostable::operator bool() const {
  return transient_;
}

}  // namespace ksp_plugin
}  // namespace principia
