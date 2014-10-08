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

inline Body::Body(GravitationalParameter const& gravitational_parameter,
                  double const j2,
                  Length const& radius,
                  R3Element<double> const& axis)
    : Body(gravitational_parameter,
           -j2 * gravitational_parameter * radius * radius,
           axis) {}

inline Body::Body(Mass const& mass,
                  double const j2,
                  Length const& radius,
                  R3Element<double> const& axis)
    : Body(mass,
           -j2 * mass * GravitationalConstant * radius * radius,
           axis) {}

inline Body::Body(GravitationalParameter const& gravitational_parameter,
                  Order2ZonalCoefficient const& j2,
                  R3Element<double> const& axis)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant),
      j2_(j2),
      axis_(axis) {}

inline Body::Body(Mass const& mass,
                  Order2ZonalCoefficient const& j2,
                  R3Element<double> const& axis)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass),
      j2_(j2),
      axis_(axis) {}

inline GravitationalParameter const& Body::gravitational_parameter() const {
  return gravitational_parameter_;
}

inline Mass const& Body::mass() const {
  return mass_;
}

inline Order2ZonalCoefficient const& Body::j2() const {
  return j2_;
}

inline R3Element<double> const& Body::axis() const {
  return axis_;
}

inline bool Body::is_massless() const {
  return mass_ == Mass();
}

inline bool Body::is_oblate() const {
  return j2_ != Order2ZonalCoefficient();
}

}  // namespace physics
}  // namespace principia
