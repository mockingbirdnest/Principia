#pragma once

#include "physics/discrete_trajectory2.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory2 {

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::begin() const {
  return (*segments_.begin())->begin();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::end() const {
  return (*segments_.rbegin())->end();
}

}  // namespace internal_discrete_trajectory2
}  // namespace physics
}  // namespace principia
