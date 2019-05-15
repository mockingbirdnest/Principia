
#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <cmath>
#include <functional>

#include "base/not_null.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_manœuvre {

using base::check_not_null;
using geometry::NormalizeOrZero;
using geometry::Rotation;
using physics::RigidMotion;
using quantities::Acceleration;
using quantities::Sqrt;
using std::placeholders::_1;

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame>::Manœuvre(
    Force const& thrust,
    Mass const& initial_mass,
    SpecificImpulse const& specific_impulse,
    Intensity const& intensity,
    Timing const& timing,
    not_null<std::shared_ptr<DynamicFrame<InertialFrame, Frame> const>> frame,
    bool is_inertially_fixed)
    : thrust_(thrust),
      initial_mass_(initial_mass),
      specific_impulse_(specific_impulse),
      intensity_(intensity),
      timing_(timing),
      frame_(std::move(frame)),
      is_inertially_fixed_(is_inertially_fixed) {
  // Fill the missing fields of |intensity_|.
  if (intensity_.Δv) {
    CHECK(!intensity_.direction && !intensity_.duration);
    intensity_.direction = NormalizeOrZero(*intensity_.Δv);
    auto const speed = intensity_.Δv->Norm();
    if (speed == Speed()) {
      // This handles the case where |thrust_| vanishes, where the usual formula
      // would yield NaN.
      intensity_.duration = Time();
    } else {
      intensity_.duration =
          initial_mass_ * specific_impulse_ *
          (1 - std::exp(-speed / specific_impulse_)) / thrust_;
    }
  } else {
    CHECK(intensity_.direction && intensity_.duration);
    // Циолко́вский's equation.
    intensity_.Δv = *intensity_.direction * specific_impulse_ *
                    std::log(initial_mass_ / final_mass());
  }

  // Fill the missing fields of |timing_|.
  if (timing_.initial_time) {
    CHECK(!timing_.time_of_half_Δv);
    timing_.time_of_half_Δv = *timing_.initial_time + time_to_half_Δv();
  } else {
    CHECK(timing_.time_of_half_Δv);
    timing_.initial_time = *timing_.time_of_half_Δv - time_to_half_Δv();
  }
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
SpecificImpulse const&
Manœuvre<InertialFrame, Frame>::specific_impulse() const {
  return specific_impulse_;
}

template<typename InertialFrame, typename Frame>
not_null<std::shared_ptr<DynamicFrame<InertialFrame, Frame> const>>
Manœuvre<InertialFrame, Frame>::frame() const {
  return frame_;
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::is_inertially_fixed() const {
  return is_inertially_fixed_;
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Intensity const&
    Manœuvre<InertialFrame, Frame>::intensity() const {
  return intensity_;
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Timing const&
    Manœuvre<InertialFrame, Frame>::timing() const {
  return timing_;
}

template<typename InertialFrame, typename Frame>
Variation<Mass> Manœuvre<InertialFrame, Frame>::mass_flow() const {
  return thrust_ / specific_impulse_;
}

template<typename InertialFrame, typename Frame>
Mass Manœuvre<InertialFrame, Frame>::final_mass() const {
  return initial_mass_ - mass_flow() * *intensity_.duration;
}

template<typename InertialFrame, typename Frame>
Time Manœuvre<InertialFrame, Frame>::time_to_half_Δv() const {
  return specific_impulse_ * initial_mass_*
         (1 - std::sqrt(final_mass() / initial_mass_)) / thrust_;
}

template<typename InertialFrame, typename Frame>
Instant Manœuvre<InertialFrame, Frame>::final_time() const {
  return *timing_.initial_time + *intensity_.duration;
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::FitsBetween(Instant const& begin,
                                                 Instant const& end) const {
  return begin < initial_time() && final_time() < end;
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::IsSingular() const {
  return !IsFinite(Δv());
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_coasting_trajectory(
    not_null<DiscreteTrajectory<InertialFrame> const*> const trajectory) {
  coasting_trajectory_ = trajectory;
}

template<typename InertialFrame, typename Frame>
Vector<double, InertialFrame>
    Manœuvre<InertialFrame, Frame>::InertialDirection() const {
  CHECK(is_inertially_fixed_);
  return FrenetFrame()(*intensity_.direction);
}

template<typename InertialFrame, typename Frame>
typename Ephemeris<InertialFrame>::IntrinsicAcceleration
Manœuvre<InertialFrame, Frame>::InertialIntrinsicAcceleration() const {
  CHECK(is_inertially_fixed_);
  return std::bind(
      &Manœuvre<InertialFrame, Frame>::ComputeIntrinsicAcceleration,
      this,
      _1,
      /*direction=*/InertialDirection());
}

template<typename InertialFrame, typename Frame>
typename Ephemeris<InertialFrame>::GeneralizedIntrinsicAcceleration
Manœuvre<InertialFrame, Frame>::FrenetIntrinsicAcceleration() const {
  CHECK(!is_inertially_fixed_);
  return [this](Instant const& t,
                DegreesOfFreedom<InertialFrame> const& degrees_of_freedom)
             -> Vector<Acceleration, InertialFrame> {
    return ComputeIntrinsicAcceleration(
        t,
        /*direction=*/ComputeFrenetFrame(t, degrees_of_freedom)(direction_));
  };
}

template<typename InertialFrame, typename Frame>
OrthogonalMap<Frenet<Frame>, InertialFrame>
    Manœuvre<InertialFrame, Frame>::FrenetFrame() const {
  CHECK_NOTNULL(coasting_trajectory_);
  typename DiscreteTrajectory<InertialFrame>::Iterator const it =
      coasting_trajectory_->Find(*timing_.initial_time);
  CHECK(it != coasting_trajectory_->End());
  return ComputeFrenetFrame(*timing_.initial_time, it.degrees_of_freedom());
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::WriteToMessage(
    not_null<serialization::Manoeuvre*> const message) const {
  thrust_.WriteToMessage(message->mutable_thrust());
  initial_mass_.WriteToMessage(message->mutable_initial_mass());
  specific_impulse_.WriteToMessage(message->mutable_specific_impulse());
  intensity_.direction->WriteToMessage(message->mutable_direction());
  intensity_.duration->WriteToMessage(message->mutable_duration());
  timing_.initial_time->WriteToMessage(message->mutable_initial_time());
  frame_->WriteToMessage(message->mutable_frame());
  message->set_is_inertially_fixed(is_inertially_fixed_);
}

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame> Manœuvre<InertialFrame, Frame>::ReadFromMessage(
    serialization::Manoeuvre const& message,
    not_null<Ephemeris<InertialFrame>*> const ephemeris) {
  Intensity intensity;
  intensity.direction =
      Vector<double, Frenet<Frame>>::ReadFromMessage(message.direction());
  intensity.duration = Time::ReadFromMessage(message.duration());
  Timing timing;
  timing.initial_time = Instant::ReadFromMessage(message.initial_time());
  return Manœuvre(
      Force::ReadFromMessage(message.thrust()),
      Mass::ReadFromMessage(message.initial_mass()),
      SpecificImpulse::ReadFromMessage(message.specific_impulse()),
      intensity,
      timing,
      DynamicFrame<InertialFrame, Frame>::ReadFromMessage(message.frame(),
                                                          ephemeris),
      message.is_inertially_fixed());
}

template<typename InertialFrame, typename Frame>
OrthogonalMap<Frenet<Frame>, InertialFrame>
Manœuvre<InertialFrame, Frame>::ComputeFrenetFrame(
    Instant const& t,
    DegreesOfFreedom<InertialFrame> const& degrees_of_freedom) const {
  RigidMotion<InertialFrame, Frame> const to_frame_at_t =
      frame_->ToThisFrameAtTime(t);
  RigidMotion<Frame, InertialFrame> const from_frame_at_t =
      to_frame_at_t.Inverse();
  return from_frame_at_t.orthogonal_map() *
         frame_->FrenetFrame(t, to_frame_at_t(degrees_of_freedom)).Forget();
}

template<typename InertialFrame, typename Frame>
Vector<Acceleration, InertialFrame>
Manœuvre<InertialFrame, Frame>::ComputeIntrinsicAcceleration(
    Instant const& t,
    Vector<double, InertialFrame> const& direction) const {
  if (t >= *timing_.initial_time && t <= final_time()) {
    return direction * thrust_ /
           (initial_mass_ - (t - *timing_.initial_time) * mass_flow());
  } else {
    return Vector<Acceleration, InertialFrame>();
  }
}

}  // namespace internal_manœuvre
}  // namespace ksp_plugin
}  // namespace principia
