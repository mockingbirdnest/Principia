#include "body.hpp"

#include <algorithm>
#include <iterator>

namespace principia {
namespace physics {

template<typename Frame>
Body<Frame>::Body(GravitationalParameter const& gravitational_parameter)
    : gravitational_parameter_(gravitational_parameter) {}

template<typename Frame>
Body<Frame>::~Body() {}

template<typename Frame>
GravitationalParameter const& Body<Frame>::gravitational_parameter() const {
  return gravitational_parameter_;
}

template<typename Frame>
bool Body<Frame>::is_massless() const {
  return gravitational_parameter_ == 0 * GravitationalParameter::SIUnit();
}

template<typename Frame>
void Body<Frame>::AppendToTrajectory(
    std::vector<Vector<Length, Frame>> const& positions,
    std::vector<Vector<Momentum, Frame>> const& momenta,
    std::vector<Time> const& times) {
  std::move(positions.begin(), positions.end(), back_inserter(positions_));
  std::move(momenta.begin(), momenta.end(), back_inserter(momenta_));
  std::move(times.begin(), times.end(), back_inserter(times_));
}

template<typename Frame>
void Body<Frame>::GetLast(Vector<Length, Frame>* position,
                          Vector<Momentum, Frame>* momentum,
                          Time* time) const {
  *position = *positions_.back();
  *momentum = *momenta_.back();
  *time = *times_.back();
}

}  // namespace physics
}  // namespace principia
