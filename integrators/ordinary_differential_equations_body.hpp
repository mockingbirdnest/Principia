
#pragma once

#include "integrators/ordinary_differential_equations.hpp"

#include <vector>

namespace principia {
namespace integrators {
namespace internal_ordinary_differential_equations {

template<typename... State>
FirstOrderOrdinaryDifferentialEquation<State...>::SystemState::SystemState(
    State const& y,
    Instant const& t)
    : y(y), time(t) {}

template<typename Position_>
SpecialSecondOrderDifferentialEquation<Position_>::SystemState::SystemState(
    std::vector<Position> const& q,
    std::vector<Velocity> const& v,
    Instant const& t)
    : time(t) {
  for (int i = 0; i < q.size(); ++i) {
    positions.emplace_back(q[i]);
    velocities.emplace_back(v[i]);
  }
}

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
  SystemState system_state;
  for (auto const& p : message.position()) {
    system_state.positions.push_back(
        DoublePrecision<Position>::ReadFromMessage(p));
  }
  for (auto const& v : message.velocity()) {
    system_state.velocities.push_back(
        DoublePrecision<Velocity>::ReadFromMessage(v));
  }
  system_state.time = DoublePrecision<Instant>::ReadFromMessage(message.time());
  return system_state;
}

}  // namespace internal_ordinary_differential_equations
}  // namespace integrators
}  // namespace principia
