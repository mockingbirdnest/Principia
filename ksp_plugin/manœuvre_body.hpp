#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <cmath>

#include "quantities/elementary_functions.hpp"

namespace principia {

using physics::RigidMotion;
using quantities::Sqrt;

namespace ksp_plugin {

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame>::Manœuvre(
    Force const& thrust,
    Mass const& initial_mass,
    SpecificImpulse const& specific_impulse,
    Vector<double, Frenet<Frame>> const& direction,
    not_null<DynamicFrame<InertialFrame, Frame> const*> frame)
    : thrust_(thrust),
      initial_mass_(initial_mass),
      specific_impulse_(specific_impulse),
      direction_(direction),
      frame_(frame) {}

template<typename InertialFrame, typename Frame>
Instant Manœuvre<InertialFrame, Frame>::initial_time() const {
  CHECK(initial_time_);
  return *initial_time_;
}

template<typename InertialFrame, typename Frame>
Instant Manœuvre<InertialFrame, Frame>::time_of_half_Δv() const {
  return initial_time() + time_to_half_Δv();
}

template<typename InertialFrame, typename Frame>
Instant Manœuvre<InertialFrame, Frame>::final_time() const {
  return initial_time() + duration();
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_initial_time(
    Instant const& initial_time) {
  initial_time_ = initial_time;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_time_of_half_Δv(
    Instant const& time_of_half_Δv) {
  set_initial_time(time_of_half_Δv - time_to_half_Δv());
}

template<typename InertialFrame, typename Frame>
Time Manœuvre<InertialFrame, Frame>::duration() const {
  CHECK(duration_);
  return *duration_;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_duration(Time const& duration) {
  duration_ = duration;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_Δv(Speed const& Δv) {
  set_duration(initial_mass_ * specific_impulse_ *
               (1 - std::exp(-Δv / specific_impulse_)) / thrust_);
}

template<typename InertialFrame, typename Frame>
Speed Manœuvre<InertialFrame, Frame>::Δv() const {
  // Циолко́вский's equation.
  return specific_impulse_ * std::log(initial_mass_ / final_mass());
}

template<typename InertialFrame, typename Frame>
Vector<double, Frenet<Frame>> const&
Manœuvre<InertialFrame, Frame>::direction() const {
  return direction_;
}

template<typename InertialFrame, typename Frame>
SpecificImpulse const&
Manœuvre<InertialFrame, Frame>::specific_impulse() const {
  return specific_impulse_;
}

template<typename InertialFrame, typename Frame>
Force const& Manœuvre<InertialFrame, Frame>::thrust() const {
  return thrust_;
}

template<typename InertialFrame, typename Frame>
Mass const& Manœuvre<InertialFrame, Frame>::initial_mass() const {
  return initial_mass_;
}

template<typename InertialFrame, typename Frame>
Variation<Mass> Manœuvre<InertialFrame, Frame>::mass_flow() const {
  return thrust_ / specific_impulse_;
}

template<typename InertialFrame, typename Frame>
Mass Manœuvre<InertialFrame, Frame>::final_mass() const {
  return initial_mass_ - mass_flow() * duration();
}

template<typename InertialFrame, typename Frame>
Time Manœuvre<InertialFrame, Frame>::time_to_half_Δv() const {
  return specific_impulse_ * initial_mass_*
         (1 - std::sqrt(final_mass() / initial_mass_)) / thrust_;
}

template<typename InertialFrame, typename Frame>
typename Ephemeris<InertialFrame>::IntrinsicAcceleration
Manœuvre<InertialFrame, Frame>::acceleration(
    DiscreteTrajectory<InertialFrame> const& coasting_trajectory) const {
  typename DiscreteTrajectory<InertialFrame>::Iterator const it =
      coasting_trajectory.Find(initial_time());
  CHECK(it != coasting_trajectory.End());
  RigidMotion<InertialFrame, Frame> const to_frame_at_initial_time =
      frame_->ToThisFrameAtTime(initial_time());
  OrthogonalMap<Frame, InertialFrame> const from_frame_at_initial_time =
      to_frame_at_initial_time.orthogonal_map().Inverse();
  Rotation<Frenet<Frame>, Frame> const from_frenet_frame = 
      frame_->FrenetFrame(initial_time(),
                          to_frame_at_initial_time(it.degrees_of_freedom()));
  Vector<double, InertialFrame> inertial_direction =
      from_frame_at_initial_time(from_frenet_frame(direction_));

  return [this, inertial_direction](
      Instant const& time) -> Vector<Acceleration, InertialFrame> {
    if (time >= initial_time() && time <= final_time()) {
      return inertial_direction * thrust_ /
             (initial_mass_ - (time - initial_time()) * mass_flow());
    } else {
      return Vector<Acceleration, InertialFrame>();
    }
  };
}

}  // namespace ksp_plugin
}  // namespace principia
