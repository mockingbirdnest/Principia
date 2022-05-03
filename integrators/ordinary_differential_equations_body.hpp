
#pragma once

#include "integrators/ordinary_differential_equations.hpp"

#include <vector>

#include "base/for_all_of.hpp"

namespace principia {
namespace integrators {
namespace internal_ordinary_differential_equations {

using base::for_all_of;

template<typename... State>
ExplicitFirstOrderOrdinaryDifferentialEquation<
    State...>::SystemState::SystemState(State const& y, Instant const& t)
    : time(t) {
  for_all_of(y, this->y).loop([](auto const& y, auto& this_y) {
    for (auto const& y_i : y) {
      this_y.emplace_back(y_i);
    }
  });
}

template<typename... State>
void ExplicitFirstOrderOrdinaryDifferentialEquation<State...>::SystemState::
    WriteToMessage(not_null<serialization::SystemState*> message) const {
  // Writing the tuple would be tricky.
  LOG(FATAL) << "NYI";
}

template<typename... State>
ExplicitFirstOrderOrdinaryDifferentialEquation<State...>::SystemState
ExplicitFirstOrderOrdinaryDifferentialEquation<State...>::SystemState::
    ReadFromMessage(serialization::SystemState const& message) {
  // Reading the tuple would be tricky.
  LOG(FATAL) << "NYI";
}

template<typename... State>
DecomposableFirstOrderDifferentialEquation<State...>::SystemState::SystemState(
    State const& y,
    Instant const& t)
    : y(y), time(t) {}

template<typename Position_>
ExplicitSecondOrderOrdinaryDifferentialEquation<
    Position_>::SystemState::SystemState(std::vector<Position> const& q,
                                         std::vector<Velocity> const& v,
                                         Instant const& t)
    : time(t) {
  for (int i = 0; i < q.size(); ++i) {
    positions.emplace_back(q[i]);
    velocities.emplace_back(v[i]);
  }
}

template<typename Position>
void ExplicitSecondOrderOrdinaryDifferentialEquation<Position>::SystemState::
    WriteToMessage(not_null<serialization::SystemState*> const message) const {
  for (auto const& position : positions) {
    position.WriteToMessage(message->add_position());
  }
  for (auto const& velocity : velocities) {
    velocity.WriteToMessage(message->add_velocity());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename Position>
typename ExplicitSecondOrderOrdinaryDifferentialEquation<Position>::SystemState
ExplicitSecondOrderOrdinaryDifferentialEquation<Position>::SystemState::
    ReadFromMessage(serialization::SystemState const& message) {
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
