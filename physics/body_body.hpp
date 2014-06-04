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
      mass_(gravitational_parameter / GravitationalConstant) {}

template<typename Frame>
Body<Frame>::Body(Mass const& mass)
    : gravitational_parameter_(mass * GravitationalConstant),
      mass_(mass) {}

template<typename Frame>
Body<Frame>::~Body() {}

template<typename Frame>
GravitationalParameter const& Body<Frame>::gravitational_parameter() const {
  return gravitational_parameter_;
}

template<typename Frame>
Mass const& Body<Frame>::mass() const {
  return mass_;
}

template<typename Frame>
bool Body<Frame>::is_massless() const {
  return mass_ == 0 * Mass::SIUnit();
}

template<typename Frame>
void Body<Frame>::AppendToTrajectory(Vector<Length, Frame> const& position,
                                     Vector<Speed, Frame> const& velocity,
                                     Time const& time) {
  positions_.push_back(position);
  velocities_.push_back(velocity);
  times_.push_back(time);
}

template<typename Frame>
std::vector<Vector<Length, Frame>> const& Body<Frame>::positions() const {
  return positions_;
}

template<typename Frame>
std::vector<Vector<Speed, Frame>> const& Body<Frame>::velocities() const {
  return velocities_;
}

template<typename Frame>
std::vector<Time> const& Body<Frame>::times() const {
  return times_;
}

}  // namespace physics
}  // namespace principia
