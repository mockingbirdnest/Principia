#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <cmath>

#include "quantities/elementary_functions.hpp"

namespace principia {

using quantities::Sqrt;

namespace ksp_plugin {

template <typename Frame>
Manœuvre<Frame>::Manœuvre(Force thrust,
                          Mass initial_mass,
                          Speed effective_exhaust_velocity,
                          Vector<double, Frame> direction)
    : thrust_(thrust),
      initial_mass_(initial_mass),
      effective_exhaust_velocity_(effective_exhaust_velocity),
      direction_(direction) {}

template<typename Frame>
Instant Manœuvre<Frame>::initial_time() const {
  CHECK(initial_time_ != nullptr);
  return *initial_time_;
}

template<typename Frame>
inline Instant Manœuvre<Frame>::time_of_half_Δv() const {
  return initial_time() + time_to_half_Δv();
}

template<typename Frame>
Instant Manœuvre<Frame>::final_time() const {
  return initial_time() + duration();
}

template<typename Frame>
void Manœuvre<Frame>::set_initial_time(Instant const & initial_time) {
  initial_time_ = std::make_unique<Instant>(initial_time);
}

template <typename Frame>
void Manœuvre<Frame>::set_time_of_half_Δv(Instant const& time_of_half_Δv) {
  set_initial_time(time_of_half_Δv - time_to_half_Δv());
}

template<typename Frame>
inline Time Manœuvre<Frame>::duration() const {
  CHECK(duration_ != nullptr);
  return *duration_;
}

template<typename Frame>
void Manœuvre<Frame>::set_duration(Time const & duration) {
  duration_ = std::make_unique<Time>(duration);
}

template<typename Frame>
void Manœuvre<Frame>::set_Δv(Speed const& Δv) {
  set_duration(initial_mass() * effective_exhaust_velocity() *
               (1 - std::exp(-Δv / effective_exhaust_velocity())) / thrust());
}

template<typename Frame>
Speed Manœuvre<Frame>::Δv() const {
  // Циолко́вский's equation.
  return effective_exhaust_velocity() * std::log(initial_mass() / final_mass());
}

template<typename Frame>
Vector<double, Frame> Manœuvre<Frame>::direction() const {
  return direction_;
}

template<typename Frame>
Speed Manœuvre<Frame>::effective_exhaust_velocity() const {
  return effective_exhaust_velocity_;
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
  return thrust_ / effective_exhaust_velocity_;
}

template<typename Frame>
Mass Manœuvre<Frame>::final_mass() const {
  return initial_mass() - mass_flow() * duration();
}

template<typename Frame>
Time Manœuvre<Frame>::time_to_half_Δv() const {
  return effective_exhaust_velocity() * initial_mass() *
         (1 - std::sqrt(final_mass() / initial_mass())) / thrust();
}

template <typename Frame>
typename Trajectory<Frame>::IntrinsicAcceleration
    Manœuvre<Frame>::acceleration() const {
  return [
    direction = this->direction(),
    initial_time = this->initial_time(),
    final_time = this->final_time(),
    thrust = this->thrust(),
    initial_mass = this->initial_mass(),
    mass_flow = this->mass_flow()
  ](Instant const& time) -> Vector<Acceleration, Frame> {
    if (time >= initial_time && time <= final_time) {
      return direction * thrust /
             (initial_mass - (time - initial_time) * mass_flow);
    } else {
      return Vector<Acceleration, Frame>();
    }
  };
}

}  // namespace ksp_plugin
}  // namespace principia
