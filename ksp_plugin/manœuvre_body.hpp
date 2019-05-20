
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
Manœuvre<InertialFrame, Frame>::Manœuvre(Mass const& initial_mass,
                                         Burn const& burn)
  : initial_mass_(initial_mass),
    construction_burn_(burn),
    burn_(burn) {
  // Fill the missing fields of |intensity|.
  auto& intensity = burn_.intensity;
  if (intensity.Δv) {
    CHECK(!intensity.direction && !intensity.duration);
    intensity.direction = NormalizeOrZero(*intensity.Δv);
    auto const speed = intensity.Δv->Norm();
    if (speed == Speed()) {
      // This handles the case where |thrust_| vanishes, where the usual formula
      // would yield NaN.
      intensity.duration = Time();
    } else {
      intensity.duration =
          initial_mass_ * specific_impulse() *
          (1 - std::exp(-speed / specific_impulse())) / thrust();
    }
  } else {
    CHECK(intensity.direction && intensity.duration);
    // Циолковский's equation.
    intensity.Δv = *intensity.direction * specific_impulse() *
                   std::log(initial_mass_ / final_mass());
  }

  // Fill the missing fields of |timing|.
  auto& timing = burn_.timing;
  if (timing.initial_time) {
    CHECK(!timing.time_of_half_Δv);
    timing.time_of_half_Δv = *timing.initial_time + time_to_half_Δv();
  } else {
    CHECK(timing.time_of_half_Δv);
    timing.initial_time = *timing.time_of_half_Δv - time_to_half_Δv();
  }
}

template<typename InertialFrame, typename Frame>
Mass const& Manœuvre<InertialFrame, Frame>::initial_mass() const {
  return initial_mass_;
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Intensity const&
    Manœuvre<InertialFrame, Frame>::intensity() const {
  return construction_burn_.intensity;
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Timing const&
    Manœuvre<InertialFrame, Frame>::timing() const {
  return construction_burn_.timing;
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Burn const&
    Manœuvre<InertialFrame, Frame>::burn() const {
  return construction_burn_;
}

template<typename InertialFrame, typename Frame>
Vector<double, Frenet<Frame>> const& Manœuvre<InertialFrame, Frame>::direction()
    const {
  return *full_intensity().direction;
}

template<typename InertialFrame, typename Frame>
Time const& Manœuvre<InertialFrame, Frame>::duration() const {
  return *full_intensity().duration;
}

template<typename InertialFrame, typename Frame>
Velocity<Frenet<Frame>> const& Manœuvre<InertialFrame, Frame>::Δv() const {
  return *full_intensity().Δv;
}

template<typename InertialFrame, typename Frame>
Instant const& Manœuvre<InertialFrame, Frame>::initial_time() const {
  return *full_timing().initial_time;
}

template<typename InertialFrame, typename Frame>
Instant const& Manœuvre<InertialFrame, Frame>::time_of_half_Δv() const {
  return *full_timing().time_of_half_Δv;
}

template<typename InertialFrame, typename Frame>
Force const& Manœuvre<InertialFrame, Frame>::thrust() const {
  return burn_.thrust;
}

template<typename InertialFrame, typename Frame>
SpecificImpulse const&
Manœuvre<InertialFrame, Frame>::specific_impulse() const {
  return burn_.specific_impulse;
}

template<typename InertialFrame, typename Frame>
not_null<std::shared_ptr<DynamicFrame<InertialFrame, Frame> const>>
Manœuvre<InertialFrame, Frame>::frame() const {
  return burn_.frame;
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::is_inertially_fixed() const {
  return burn_.is_inertially_fixed;
}

template<typename InertialFrame, typename Frame>
Variation<Mass> Manœuvre<InertialFrame, Frame>::mass_flow() const {
  return thrust() / specific_impulse();
}

template<typename InertialFrame, typename Frame>
Mass Manœuvre<InertialFrame, Frame>::final_mass() const {
  return initial_mass_ - mass_flow() * duration();
}

template<typename InertialFrame, typename Frame>
Time Manœuvre<InertialFrame, Frame>::time_to_half_Δv() const {
  return specific_impulse() * initial_mass_*
         (1 - std::sqrt(final_mass() / initial_mass_)) / thrust();
}

template<typename InertialFrame, typename Frame>
Instant Manœuvre<InertialFrame, Frame>::final_time() const {
  return initial_time() + duration();
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::FitsBetween(Instant const& begin,
                                                 Instant const& end) const {
  return begin < initial_time() && final_time() < end;
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::IsSingular() const {
  return !IsFinite(Δv().Norm²());
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_coasting_trajectory(
    not_null<DiscreteTrajectory<InertialFrame> const*> const trajectory) {
  coasting_trajectory_ = trajectory;
}

template<typename InertialFrame, typename Frame>
Vector<double, InertialFrame>
    Manœuvre<InertialFrame, Frame>::InertialDirection() const {
  CHECK(is_inertially_fixed());
  return FrenetFrame()(direction());
}

template<typename InertialFrame, typename Frame>
typename Ephemeris<InertialFrame>::IntrinsicAcceleration
Manœuvre<InertialFrame, Frame>::InertialIntrinsicAcceleration() const {
  CHECK(is_inertially_fixed());
  return std::bind(
      &Manœuvre<InertialFrame, Frame>::ComputeIntrinsicAcceleration,
      this,
      _1,
      /*direction=*/InertialDirection());
}

template<typename InertialFrame, typename Frame>
typename Ephemeris<InertialFrame>::GeneralizedIntrinsicAcceleration
Manœuvre<InertialFrame, Frame>::FrenetIntrinsicAcceleration() const {
  CHECK(!is_inertially_fixed());
  return [this](Instant const& t,
                DegreesOfFreedom<InertialFrame> const& degrees_of_freedom)
             -> Vector<Acceleration, InertialFrame> {
    return ComputeIntrinsicAcceleration(
        t,
        /*direction=*/ComputeFrenetFrame(t, degrees_of_freedom)(direction()));
  };
}

template<typename InertialFrame, typename Frame>
OrthogonalMap<Frenet<Frame>, InertialFrame>
    Manœuvre<InertialFrame, Frame>::FrenetFrame() const {
  if (coasting_trajectory_==nullptr) {
    LOG(FATAL)<<"oops";
  }
  CHECK_NOTNULL(coasting_trajectory_);
  typename DiscreteTrajectory<InertialFrame>::Iterator const it =
      coasting_trajectory_->Find(initial_time());
  CHECK(it != coasting_trajectory_->End());
  return ComputeFrenetFrame(initial_time(), it.degrees_of_freedom());
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::WriteToMessage(
    not_null<serialization::Manoeuvre*> const message) const {
  thrust().WriteToMessage(message->mutable_thrust());
  initial_mass_.WriteToMessage(message->mutable_initial_mass());
  specific_impulse().WriteToMessage(message->mutable_specific_impulse());
  direction().WriteToMessage(message->mutable_direction());
  duration().WriteToMessage(message->mutable_duration());
  initial_time().WriteToMessage(message->mutable_initial_time());
  frame()->WriteToMessage(message->mutable_frame());
  message->set_is_inertially_fixed(is_inertially_fixed());
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
  Burn const burn{intensity,
                  timing,
                  Force::ReadFromMessage(message.thrust()),
                  SpecificImpulse::ReadFromMessage(message.specific_impulse()),
                  DynamicFrame<InertialFrame, Frame>::ReadFromMessage(
                      message.frame(), ephemeris),
                  message.is_inertially_fixed()};
  return Manœuvre(Mass::ReadFromMessage(message.initial_mass()), burn);
}

template<typename InertialFrame, typename Frame>
OrthogonalMap<Frenet<Frame>, InertialFrame>
Manœuvre<InertialFrame, Frame>::ComputeFrenetFrame(
    Instant const& t,
    DegreesOfFreedom<InertialFrame> const& degrees_of_freedom) const {
  RigidMotion<InertialFrame, Frame> const to_frame_at_t =
      frame()->ToThisFrameAtTime(t);
  RigidMotion<Frame, InertialFrame> const from_frame_at_t =
      to_frame_at_t.Inverse();
  return from_frame_at_t.orthogonal_map() *
         frame()->FrenetFrame(t, to_frame_at_t(degrees_of_freedom)).Forget();
}

template<typename InertialFrame, typename Frame>
Vector<Acceleration, InertialFrame>
Manœuvre<InertialFrame, Frame>::ComputeIntrinsicAcceleration(
    Instant const& t,
    Vector<double, InertialFrame> const& direction) const {
  if (t >= initial_time() && t <= final_time()) {
    return direction * thrust() /
           (initial_mass_ - (t - initial_time()) * mass_flow());
  } else {
    return Vector<Acceleration, InertialFrame>();
  }
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Intensity const&
Manœuvre<InertialFrame, Frame>::full_intensity() const {
  return burn_.intensity;
}

template<typename InertialFrame, typename Frame>
typename Manœuvre<InertialFrame, Frame>::Timing const&
Manœuvre<InertialFrame, Frame>::full_timing() const {
  return burn_.timing;
}

}  // namespace internal_manœuvre
}  // namespace ksp_plugin
}  // namespace principia
