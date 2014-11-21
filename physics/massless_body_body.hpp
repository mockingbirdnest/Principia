#pragma once

#include "physics/massless_body.hpp"

namespace principia {
namespace physics {

inline bool MasslessBody::is_massless() const {
  return true;
}

inline bool MasslessBody::is_oblate() const {
  return false;
}

}  // namespace physics
}  // namespace principia
