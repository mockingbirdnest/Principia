#pragma once

#include "physics/massive_body.hpp"

#include "glog/logging.h"
#include "quantities/constants.hpp"

using principia::constants::GravitationalConstant;

namespace principia {
namespace physics {

MassiveBody::MassiveBody(GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant) {
  CHECK_NE(gravitational_parameter, GravitationalParameter())
      << "Massive body cannot have zero gravitational parameter";
}

MassiveBody::MassiveBody(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {
  CHECK_NE(mass, Mass())
      << "Massive body cannot have zero mass";
}

template<typename Frame>
bool MassiveBody::is_massless() const {
  return false;
}

template<typename Frame>
bool MassiveBody::is_oblate() const {
  return false;
}

}  // namespace physics
}  // namespace principia
