#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <cmath>

#include "quantities/elementary_functions.hpp"

namespace principia {

using quantities::Sqrt;

namespace ksp_plugin {

template<typename Frame>
Manœuvre<Frame>::Manœuvre(Force const& thrust,
                          Mass const& initial_mass,
                          SpecificImpulse const& specific_impulse,
                          Vector<double, Frame> const& direction)
    : thrust_(thrust),
      initial_mass_(initial_mass),
      specific_impulse_(specific_impulse),
      direction_(direction) {}

template<typename Frame>
Instant Manœuvre<Frame>::initial_time() const {
  CHECK(initial_time_);
  return *initial_time_;
}

template<typename Frame>
Instant Manœuvre<Frame>::time_of_half_Δv() const {
  return initial_time() + time_to_half_Δv();
}

template<typename Frame>
Instant Manœuvre<Frame>::final_time() const {
  return initial_time() + duration();
}

template<typename Frame>
void Manœuvre<Frame>::set_initial_time(Instant const& initial_time) {
  initial_time_ = initial_time;
}

template <typename Frame>
void Manœuvre<Frame>::set_time_of_half_Δv(Instant const& time_of_half_Δv) {
  set_initial_time(time_of_half_Δv - time_to_half_Δv());
}

template<typename Frame>
Time Manœuvre<Frame>::duration() const {
  CHECK(duration_);
  return *duration_;
}

template<typename Frame>
void Manœuvre<Frame>::set_duration(Time const& duration) {
  duration_ = duration;
}

template<typename Frame>
void Manœuvre<Frame>::set_Δv(Speed const& Δv) {
  set_duration(initial_mass() * specific_impulse() *
               (1 - std::exp(-Δv / specific_impulse())) / thrust());
}

template<typename Frame>
Speed Manœuvre<Frame>::Δv() const {
  // Циолко́вский's equation.
  return specific_impulse() * std::log(initial_mass() / final_mass());
}

template<typename Frame>
Vector<double, Frame> Manœuvre<Frame>::direction() const {
  return direction_;
}

template<typename Frame>
SpecificImpulse Manœuvre<Frame>::specific_impulse() const {
  return specific_impulse_;
}

template<typename Frame>
Force Manœuvre<Frame>::thrust() const {
  return thrust_;
}

template<typename Frame>
Mass Manœuvre<Frame>::initial_mass() const {
  return initial_mass_;
}

template<typename Frame>
Variation<Mass> Manœuvre<Frame>::mass_flow() const {
  return thrust() / specific_impulse();
}

template<typename Frame>
Mass Manœuvre<Frame>::final_mass() const {
  return initial_mass() - mass_flow() * duration();
}

template<typename Frame>
Time Manœuvre<Frame>::time_to_half_Δv() const {
  return specific_impulse() * initial_mass() *
         (1 - std::sqrt(final_mass() / initial_mass())) / thrust();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::IntrinsicAcceleration
    Manœuvre<Frame>::acceleration() const {
  return [this](Instant const& time) -> Vector<Acceleration, Frame> {
    if (time >= initial_time() && time <= final_time()) {
      return direction() * thrust() /
             (initial_mass() - (time - initial_time()) * mass_flow());
    } else {
      return Vector<Acceleration, Frame>();
    }
  };
}

}  // namespace ksp_plugin
}  // namespace principia
