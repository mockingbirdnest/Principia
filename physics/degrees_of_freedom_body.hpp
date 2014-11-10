#pragma once

#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {

template<typename Frame>
DegreesOfFreedom<Frame>::DegreesOfFreedom(Position<Frame> const& position,
                                          Velocity<Frame> const& velocity)
    : position(position),
      velocity(velocity) {}

template<typename Frame>
bool operator==(DegreesOfFreedom<Frame> const& left,
                DegreesOfFreedom<Frame> const& right) {
  return left.position == right.position &&
         left.velocity == right.velocity;
}

template<typename Frame, typename Weight>
DegreesOfFreedom<Frame> Barycentre(
    std::vector<DegreesOfFreedom<Frame>> const& degrees_of_freedom,
    std::vector<Weight> const& weights) {
  // This weird expression is 0, but it defines the type of the variable.
  auto positions_weighted_sum =
      (degrees_of_freedom[0].position -
       degrees_of_freedom[0].position).coordinates() * weights[0];
  auto velocities_weighted_sum =
      degrees_of_freedom[0].velocity.coordinates() * weights[0];
  Weight weight = weights[0];
  for (size_t i = 1; i < degrees_of_freedom.size(); ++i) {
    positions_weighted_sum +=
        (degrees_of_freedom[i].position -
         degrees_of_freedom[0].position).coordinates() * weights[i];
    velocities_weighted_sum +=
        degrees_of_freedom[i].velocity.coordinates() * weights[i];
    weight += weights[i];
  }
  return {degrees_of_freedom[0].position +
              Displacement<Frame>(positions_weighted_sum / weight),
          Velocity<Frame>(velocities_weighted_sum / weight)};
}

}  // namespace physics
}  // namespace principia
