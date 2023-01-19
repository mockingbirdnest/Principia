#pragma once

#include <type_traits>

#include "physics/integration_parameters.hpp"

namespace principia {
namespace physics {
namespace internal_integration_parameters {

template<typename ODE>
template<typename E, std::enable_if_t<E::order == 1, std::nullptr_t>>
AdaptiveStepParameters<ODE>::AdaptiveStepParameters(
    AdaptiveStepSizeIntegrator<E> const& integrator,
    std::int64_t const max_steps,
    Length const& length_integration_tolerance)
    : integrator_(&integrator),
      max_steps_(max_steps),
      length_integration_tolerance_(length_integration_tolerance),
      speed_integration_tolerance_(std::nullopt) {
  CHECK_LT(0, max_steps_);
  CHECK_LT(Length(), length_integration_tolerance_);
}

template<typename ODE>
template<typename E, std::enable_if_t<E::order == 2, std::nullptr_t>>
AdaptiveStepParameters<ODE>::AdaptiveStepParameters(
    AdaptiveStepSizeIntegrator<E> const& integrator,
    std::int64_t const max_steps,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance)
    : integrator_(&integrator),
      max_steps_(max_steps),
      length_integration_tolerance_(length_integration_tolerance),
      speed_integration_tolerance_(speed_integration_tolerance) {
  CHECK_LT(0, max_steps_);
  CHECK_LT(Length(), length_integration_tolerance_);
  CHECK_LT(Speed(), speed_integration_tolerance_.value());
}

template<typename ODE>
AdaptiveStepSizeIntegrator<ODE> const&
AdaptiveStepParameters<ODE>::integrator() const {
  return *integrator_;
}

template<typename ODE>
std::int64_t AdaptiveStepParameters<ODE>::max_steps() const {
  return max_steps_;
}

template<typename ODE>
Length AdaptiveStepParameters<ODE>::length_integration_tolerance() const {
  return length_integration_tolerance_;
}

template<typename ODE>
template<typename E, std::enable_if_t<E::order == 2, std::nullptr_t>>
Speed AdaptiveStepParameters<ODE>::speed_integration_tolerance() const {
  return speed_integration_tolerance_.value();
}

template<typename ODE>
void AdaptiveStepParameters<ODE>::set_max_steps(std::int64_t const max_steps) {
  CHECK_LT(0, max_steps);
  max_steps_ = max_steps;
}

template<typename ODE>
void AdaptiveStepParameters<ODE>::set_length_integration_tolerance(
    Length const& length_integration_tolerance) {
  length_integration_tolerance_ = length_integration_tolerance;
}

template<typename ODE>
template<typename E, std::enable_if_t<E::order == 2, std::nullptr_t>>
void AdaptiveStepParameters<ODE>::set_speed_integration_tolerance(
    Speed const& speed_integration_tolerance) {
  speed_integration_tolerance_ = speed_integration_tolerance;
}

template<typename ODE>
void AdaptiveStepParameters<ODE>::WriteToMessage(
    not_null<serialization::AdaptiveStepParameters*> const message) const {
  integrator_->WriteToMessage(message->mutable_integrator());
  message->set_max_steps(max_steps_);
  length_integration_tolerance_.WriteToMessage(
      message->mutable_length_integration_tolerance());
  if (speed_integration_tolerance_.has_value()) {
    speed_integration_tolerance_.value().WriteToMessage(
        message->mutable_speed_integration_tolerance());
  }
}

template<typename ODE>
AdaptiveStepParameters<ODE> AdaptiveStepParameters<ODE>::ReadFromMessage(
    serialization::AdaptiveStepParameters const& message) {
  static_assert(ODE::order == 1 || ODE::order == 2);
  if constexpr (ODE::order == 1) {
    CHECK(!message.has_speed_integration_tolerance())
        << "Unable to read from message " << message.DebugString();
    return AdaptiveStepParameters(
        AdaptiveStepSizeIntegrator<ODE>::ReadFromMessage(message.integrator()),
        message.max_steps(),
        Length::ReadFromMessage(message.length_integration_tolerance()));
  } else if constexpr (ODE::order == 2) {
    CHECK(message.has_speed_integration_tolerance())
        << "Unable to read from message " << message.DebugString();
    return AdaptiveStepParameters(
        AdaptiveStepSizeIntegrator<ODE>::ReadFromMessage(message.integrator()),
        message.max_steps(),
        Length::ReadFromMessage(message.length_integration_tolerance()),
        Speed::ReadFromMessage(message.speed_integration_tolerance()));
  }
}

template<typename ODE>
FixedStepParameters<ODE>::FixedStepParameters(
    FixedStepSizeIntegrator<ODE> const& integrator,
    Time const& step)
    : integrator_(&integrator), step_(step) {
  CHECK_LT(Time(), step);
}

template<typename ODE>
FixedStepSizeIntegrator<ODE> const&
FixedStepParameters<ODE>::integrator() const {
  return *integrator_;
}

template<typename ODE>
Time const& FixedStepParameters<ODE>::step() const {
  return step_;
}

template<typename ODE>
void FixedStepParameters<ODE>::WriteToMessage(
    not_null<serialization::FixedStepParameters*> const message) const {
  integrator_->WriteToMessage(message->mutable_integrator());
  step_.WriteToMessage(message->mutable_step());
}

template<typename ODE>
FixedStepParameters<ODE> FixedStepParameters<ODE>::ReadFromMessage(
    serialization::FixedStepParameters const& message) {
  return FixedStepParameters(
      FixedStepSizeIntegrator<ODE>::ReadFromMessage(message.integrator()),
      Time::ReadFromMessage(message.step()));
}

}  // namespace internal_integration_parameters
}  // namespace physics
}  // namespace principia
