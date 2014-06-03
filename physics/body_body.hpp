#include "body.hpp"

#include <algorithm>
#include <iterator>

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
void Body<Frame>::AppendToTrajectory(
    std::vector<Vector<Length, Frame>> const& positions,
    std::vector<Vector<Speed, Frame>> const& velocities,
    std::vector<Time> const& times) {
  std::move(positions.begin(), positions.end(), back_inserter(positions_));
  std::move(velocities.begin(), velocities.end(), back_inserter(velocities_));
  std::move(times.begin(), times.end(), back_inserter(times_));
}

template<typename Frame>
void Body<Frame>::GetTrajectory(std::vector<Vector<Length, Frame>>* positions,
                                std::vector<Vector<Speed, Frame>>* velocities,
                                std::vector<Time>* times) {
  // TODO(phl): Avoid copies here.
  positions->assign(positions_.begin(), positions_.end());
  velocities->assign(velocities_.begin(), velocities_.end());
  times->assign(times_.begin(), times_.end());
}

template<typename Frame>
void Body<Frame>::GetLast(Vector<Length, Frame>* position,
                          Vector<Speed, Frame>* velocity,
                          Time* time) const {
  *position = positions_.back();
  *velocity = velocities_.back();
  *time = times_.back();
}

}  // namespace physics
}  // namespace principia
