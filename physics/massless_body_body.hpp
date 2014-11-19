#include "body.hpp"

#include <algorithm>
#include <vector>

#include "quantities/constants.hpp"

using principia::constants::GravitationalConstant;

namespace principia {
namespace physics {

template<typename Frame>
Body<Frame>::Body(GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant),
      axis_({0, 0, 0}) {}

template<typename Frame>
Body<Frame>::Body(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass),
      axis_({0, 0, 0}) {}

template<typename Frame>
template<typename F>
Body<Frame>::Body(
    GravitationalParameter const& gravitational_parameter,
    double const j2,
    Length const& radius,
    std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis)
    : Body(gravitational_parameter,
           -j2 * gravitational_parameter * radius * radius,
           axis) {}

template<typename Frame>
template<typename F>
Body<Frame>::Body(
    Mass const& mass,
    double const j2,
    Length const& radius,
    std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis)
    : Body(mass,
           -j2 * mass * GravitationalConstant * radius * radius,
           axis) {}

template<typename Frame>
template<typename F>
Body<Frame>::Body(
    GravitationalParameter const& gravitational_parameter,
    Order2ZonalCoefficient const& j2,
    std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis)
    : gravitational_parameter_(gravitational_parameter),
      mass_(gravitational_parameter / GravitationalConstant),
      j2_(j2),
      axis_(axis) {}

template<typename Frame>
template<typename F>
Body<Frame>::Body(
    Mass const& mass,
    Order2ZonalCoefficient const& j2,
    std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass),
      j2_(j2),
      axis_(axis) {}

template<typename Frame>
GravitationalParameter const& Body<Frame>::gravitational_parameter() const {
  return gravitational_parameter_;
}

template<typename Frame>
Mass const& Body<Frame>::mass() const {
  return mass_;
}

template<typename Frame>
Order2ZonalCoefficient const& Body<Frame>::j2() const {
  return j2_;
}

template<typename Frame>
Vector<double, Frame> const& Body<Frame>::axis() const {
  return axis_;
}

template<typename Frame>
bool Body<Frame>::is_massless() const {
  return mass_ == Mass();
}

template<typename Frame>
bool Body<Frame>::is_oblate() const {
  return j2_ != Order2ZonalCoefficient();
}

}  // namespace physics
}  // namespace principia
