#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <cmath>

#include "quantities/elementary_functions.hpp"

namespace principia {

using quantities::Sqrt;

namespace ksp_plugin {

template<typename Frame>
Instant Manœuvre<Frame>::start_time() const {
  return start_time_;
}

template<typename Frame>
inline Instant Manœuvre<Frame>::half_?v() const {
  double sqrt_mass_ratio = Sqrt(mass_ratio());
  return start_time_ +
         exhaust_velocity_ * initial_mass_ * (sqrt_mass_ratio - 1) /
             (thrust_ * sqrt_mass_ratio);
}

template<typename Frame>
Instant Manœuvre<Frame>::end_time() const {
  return end_time_;
}

template<typename Frame>
inline Time Manœuvre<Frame>::duration() const {
  return end_time_ - start_time_;
}

template<typename Frame>
Speed Manœuvre<Frame>::?v() const {
  // Циолко́вский's equation.
  return exhaust_velocity_ * std::log(mass_ratio());
}

template<typename Frame>
Vector<double, Frame> Manœuvre<Frame>::direction() const {
  return direction_;
}

template<typename Frame>
Speed Manœuvre<Frame>::exhaust_velocity() const {
  return exhaust_velocity_;
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
  return thrust_ / exhaust_velocity_;
}

template<typename Frame>
Mass Manœuvre<Frame>::final_mass() const {
  return initial_mass() - mass_flow() * duration();
}

template<typename Frame>
inline double Manœuvre<Frame>::mass_ratio() const {
  return initial_mass() / final_mass();
}

template <typename Frame>
typename Trajectory<Frame>::IntrinsicAcceleration
    Manœuvre<Frame>::acceleration() const {
  return [
    direction = direction_,
    start_time = start_time_,
    end_time = end_time_,
    thrust = thrust_,
    initial_mass = initial_mass_,
    mass_flow = this->mass_flow()
  ](Instant const& time) {
    if (time > start_time && time < end_time) {
      return direction * thrust /
             (initial_mass - (time - start_time) * mass_flow);
    } else {
      return Vector<Acceleration, Frame>();
    }
  };
}

}  // namespace ksp_plugin
}  // namespace principia
