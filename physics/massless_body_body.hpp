#pragma once

#include "physics/massless_body.hpp"

namespace principia {
namespace physics {

bool MasslessBody::is_massless() const {
  return true;
}

bool MasslessBody::is_oblate() const {
  return false;
}

}  // namespace physics
}  // namespace principia
