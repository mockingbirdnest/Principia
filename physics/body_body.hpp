#include "body.hpp"

#include <algorithm>
#include <vector>

#include "quantities/constants.hpp"

// TODO(phl): Polluting the root namespace!!!
using principia::constants::GravitationalConstant;

namespace principia {
namespace physics {

Body::Body(GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant) {}

Body::Body(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {}

GravitationalParameter const& Body::gravitational_parameter() const {
  return gravitational_parameter_;
}

Mass const& Body::mass() const {
  return mass_;
}

bool Body::is_massless() const {
  return mass_ == Mass();
}

}  // namespace physics
}  // namespace principia
