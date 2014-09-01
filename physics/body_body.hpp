#include "body.hpp"

#include <algorithm>
#include <vector>

#include "quantities/constants.hpp"

using principia::constants::GravitationalConstant;

namespace principia {
namespace physics {

inline Body::Body(GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant) {}

inline Body::Body(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {}

inline GravitationalParameter const& Body::gravitational_parameter() const {
  return gravitational_parameter_;
}

inline Mass const& Body::mass() const {
  return mass_;
}

inline bool Body::is_massless() const {
  return mass_ == Mass();
}

}  // namespace physics
}  // namespace principia
