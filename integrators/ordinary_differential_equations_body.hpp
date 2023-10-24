#pragma once

#include "integrators/ordinary_differential_equations.hpp"

#include <vector>

#include "base/for_all_of.hpp"  // 🧙 For for_all_of.

namespace principia {
namespace integrators {
namespace _ordinary_differential_equations {
namespace internal {

template<typename IndependentVariable_, typename... DependentVariable>
ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable_,
                                               DependentVariable...>::
State::State(IndependentVariable const& s, DependentVariables const& y)
    : s(s) {
  using principia::base::_for_all_of::for_all_of;
  for_all_of(y, this->y).loop([](auto const& y, auto& this_y) {
    this_y = DoublePrecision(y);
  });
}

template<typename IndependentVariable_, typename... DependentVariable>
void ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable_,
                                                    DependentVariable...>::
State::WriteToMessage(not_null<serialization::State*> message) const {
  // Writing the tuple would be tricky.
  LOG(FATAL) << "NYI";
}

template<typename IndependentVariable_, typename... DependentVariable>
typename ExplicitFirstOrderOrdinaryDifferentialEquation<
    IndependentVariable_,
    DependentVariable...>::State
ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable_,
                                               DependentVariable...>::
State::ReadFromMessage(serialization::State const& message) {
  // Reading the tuple would be tricky.
  LOG(FATAL) << "NYI";
}

template<typename... DependentVariable>
DecomposableFirstOrderDifferentialEquation<DependentVariable...>::
State::State(IndependentVariable const& t, DependentVariables const& y)
    : time(t), y(y) {}

template<typename DependentVariable>
ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
State::State(IndependentVariable const& t,
             DependentVariables const& q,
             DependentVariableDerivatives const& v)
    : time(t) {
  for (int i = 0; i < q.size(); ++i) {
    positions.emplace_back(q[i]);
    velocities.emplace_back(v[i]);
  }
}

template<typename DependentVariable>
void ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
State::WriteToMessage(not_null<serialization::State*> const message) const {
  for (auto const& position : positions) {
    position.WriteToMessage(message->add_position());
  }
  for (auto const& velocity : velocities) {
    velocity.WriteToMessage(message->add_velocity());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename DependentVariable>
typename ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
    State
ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
State::ReadFromMessage(serialization::State const& message) {
  State state;
  for (auto const& p : message.position()) {
    state.positions.push_back(
        DoublePrecision<DependentVariable>::ReadFromMessage(p));
  }
  for (auto const& v : message.velocity()) {
    state.velocities.push_back(
        DoublePrecision<DependentVariableDerivative>::ReadFromMessage(v));
  }
  state.time =
      DoublePrecision<IndependentVariable>::ReadFromMessage(message.time());
  return state;
}

}  // namespace internal
}  // namespace _ordinary_differential_equations
}  // namespace integrators
}  // namespace principia
