#pragma once

#include "integrators/ordinary_differential_equations.hpp"

namespace principia {
namespace integrators {

template<typename Position>
void
SpecialSecondOrderDifferentialEquation<Position>::SystemState::WriteToMessage(
        not_null<serialization::SystemState*> const message) const {
  for (auto const& position : positions) {
    position.WriteToMessage(message->add_position());
  }
  for (auto const& velocity : velocities) {
    velocity.WriteToMessage(message->add_velocity());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename Position>
typename SpecialSecondOrderDifferentialEquation<Position>::SystemState
SpecialSecondOrderDifferentialEquation<Position>::SystemState::ReadFromMessage(
        serialization::SystemState const& message) {
  for (auto const p : message.position()) {
    positions.push_back(DoublePrecision<Position>::ReadFromMessage(p));
  }
  for (auto const v : message.velocity()) {
    velocities.push_back(DoublePrecision<Velocity>::ReadFromMessage(v));
  }
  time = Instant::ReadFromMessage(message.time());
}

}  // namespace integrators
}  // namespace principia
