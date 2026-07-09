#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <functional>
#include <memory>

#include "physics/discrete_trajectory.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace ksp_plugin {
namespace _manœuvre {
namespace internal {

using std::placeholders::_1;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_rigid_motion;

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame>::Intensity::Intensity(
    R3Element<Speed> const& Δv_cartesian_coordinates)
    : Δv_coordinates_(Δv_cartesian_coordinates),
      Δv_(Δv_cartesian_coordinates),
      direction_(NormalizeOrZero(Δv_)) {}

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame>::Intensity::Intensity(
    EvenPermutation permutation,
    SphericalCoordinates<Speed> const& Δv_spherical_coordinates)
    : Δv_coordinates_(SphericalIntensity{
          Permutation<PermutedFrenet<Frame>, Frenet<Frame>>(permutation),
          Δv_spherical_coordinates}),
      Δv_(std::get<SphericalIntensity>(Δv_coordinates_)
              .permutation(Velocity<PermutedFrenet<Frame>>(
                  Δv_spherical_coordinates.ToCartesian()))),
      direction_(NormalizeOrZero(Δv_)) {}

template<typename InertialFrame, typename Frame>
Vector<double, Frenet<Frame>> const&
Manœuvre<InertialFrame, Frame>::Intensity::direction() const {
  return direction_;
}

template<typename InertialFrame, typename Frame>
Velocity<Frenet<Frame>> const&
Manœuvre<InertialFrame, Frame>::Intensity::Δv() const {
  return Δv_;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::Intensity::set_Δv(
    Velocity<Frenet<Frame>> const& Δv) {
  if (std::holds_alternative<R3Element<Speed>>(Δv_coordinates_)) {
    Δv_coordinates_ = Δv.coordinates();
  } else {
    auto& spherical_intensity = std::get<SphericalIntensity>(Δv_coordinates_);
    Velocity<PermutedFrenet<Frame>> const permuted_Δv =
        spherical_intensity.permutation.Inverse()(Δv);
    spherical_intensity.Δv_spherical_coordinates =
        permuted_Δv.coordinates().ToSpherical();
    // The `permutation` is unchanged.
  }
  Δv_ = Δv;
  direction_ = NormalizeOrZero(Δv_);
}

template<typename InertialFrame, typename Frame>
bool
Manœuvre<InertialFrame, Frame>::Intensity::has_spherical_coordinates() const {
  return std::holds_alternative<SphericalIntensity>(Δv_coordinates_);
}

template<typename InertialFrame, typename Frame>
R3Element<Speed> const& Manœuvre<
    InertialFrame, Frame>::Intensity::Δv_cartesian_coordinates() const {
  CHECK(has_spherical_coordinates());
  return std::get<R3Element<Speed>>(Δv_coordinates_);
}

template<typename InertialFrame, typename Frame>
Permutation<PermutedFrenet<Frame>, Frenet<Frame>> const&
Manœuvre<InertialFrame, Frame>::Intensity::permutation() const {
  CHECK(has_spherical_coordinates());
  return std::get<SphericalIntensity>(Δv_coordinates_).permutation;
}

template<typename InertialFrame, typename Frame>
SphericalCoordinates<Speed> const&
Manœuvre<InertialFrame, Frame>::Intensity::Δv_spherical_coordinates() const {
  CHECK(has_spherical_coordinates());
  return std::get<SphericalIntensity>(Δv_coordinates_).Δv_spherical_coordinates;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::Intensity::WriteToMessage(
    not_null<serialization::Intensity*> const message) const {
  if (std::holds_alternative<R3Element<Speed>>(Δv_coordinates_)) {
    Δv_cartesian_coordinates().WriteToMessage(message->mutable_cartesian());
  } else {
    not_null<serialization::Intensity::Spherical*> const spherical_message =
        message->mutable_spherical();
    Δv_spherical_coordinates().WriteToMessage(
        spherical_message->mutable_coordinates());
    permutation().WriteToMessage(spherical_message->mutable_permutation());
    auto const& spherical_intensity =
        std::get<SphericalIntensity>(Δv_coordinates_);
  }
}

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame>::Intensity
Manœuvre<InertialFrame, Frame>::Intensity::ReadFromMessage(
    serialization::Intensity const& message) {
  switch (message.intensity_case()) {
    case serialization::Intensity::IntensityCase::kCartesian:
      return Intensity(R3Element<Speed>::ReadFromMessage(message.cartesian()));
    case serialization::Intensity::IntensityCase::kSpherical:
      // FIXME The proto's permutation encodes the coordinate permutation
      // integer;
      //  cast it to our EvenPermutation enum.  Spherical intensities are
      //  defined using an even permutation of (T, N, B).
      {
        auto const& spherical_message = message.spherical();
        auto const permutation_int =
            spherical_message.permutation().coordinate_permutation();
        EvenPermutation const permutation =
            static_cast<EvenPermutation>(permutation_int);
        SphericalCoordinates<Speed> const coordinates =
            SphericalCoordinates<Speed>::ReadFromMessage(
                spherical_message.coordinates());
        return Intensity(permutation, coordinates);
      }
    case serialization::Intensity::IntensityCase::INTENSITY_NOT_SET:
      LOG(FATAL) << "Missing intensity: " << message;
  };
#if PRINCIPIA_COMPILER_MSVC && \
    (_MSC_FULL_VER == 194'435'222 || _MSC_FULL_VER == 194'435'224)
  std::abort();
#endif
}

template<typename InertialFrame, typename Frame>
Manœuvre<InertialFrame, Frame>::Manœuvre(Mass const& initial_mass,
                                         Burn const& burn)
  : initial_mass_(initial_mass),
    construction_burn_(burn),
    burn_(burn) {
  auto const speed = burn.intensity.Δv().Norm();
  if (speed == Speed()) {
    // This handles the case where `thrust_` vanishes, where the usual formula
    // would yield NaN.
    duration_ = Time();
  } else {
    duration_ = initial_mass_ * specific_impulse() *
                (1 - std::exp(-speed / specific_impulse())) / thrust();
  }

  // Fill the missing fields of `timing`.
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
  return burn_.intensity.direction();
}

template<typename InertialFrame, typename Frame>
Time const& Manœuvre<InertialFrame, Frame>::duration() const {
  return duration_;
}

template<typename InertialFrame, typename Frame>
Velocity<Frenet<Frame>> const& Manœuvre<InertialFrame, Frame>::Δv() const {
  return burn_.intensity.Δv();
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
not_null<std::shared_ptr<RigidReferenceFrame<InertialFrame, Frame> const>>
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
bool Manœuvre<InertialFrame, Frame>::IsSingular() const {
  return !IsFinite(Δv().Norm²());
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::IsAfter(Instant const& time) const {
  return time < initial_time();
}

template<typename InertialFrame, typename Frame>
bool Manœuvre<InertialFrame, Frame>::FitsBetween(Instant const& begin,
                                                 Instant const& end) const {
  return begin < initial_time() && final_time() < end;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::clear_coasting_trajectory() {
  initial_degrees_of_freedom_ = std::nullopt;
}

template<typename InertialFrame, typename Frame>
void Manœuvre<InertialFrame, Frame>::set_coasting_trajectory(
    DiscreteTrajectorySegmentIterator<InertialFrame> const trajectory) {
  typename DiscreteTrajectory<InertialFrame>::iterator const it =
      trajectory->find(initial_time());
  CHECK(it != trajectory->end());
  initial_degrees_of_freedom_ = it->degrees_of_freedom;
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
  CHECK(initial_degrees_of_freedom_.has_value());
  return ComputeFrenetFrame(initial_time(),
                            initial_degrees_of_freedom_.value());
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
  bool const is_pre_leray = message.has_direction() || message.has_duration();
  LOG_IF(WARNING, is_pre_leray) << "Reading pre-Leray Manœuvre";

  Timing timing;
  timing.initial_time = Instant::ReadFromMessage(message.initial_time());
  Force const thrust = Force::ReadFromMessage(message.thrust());
  SpecificImpulse const specific_impulse =
      SpecificImpulse::ReadFromMessage(message.specific_impulse());
  Mass const initial_mass = Mass::ReadFromMessage(message.initial_mass());

  std::optional<Intensity> intensity;
  if (is_pre_leray) {
    auto const direction =
        Vector<double, Frenet<Frame>>::ReadFromMessage(message.direction());
    auto const duration = Time::ReadFromMessage(message.duration());
    auto const speed = ComputeЦиолковскийSpeed(
        initial_mass, duration, thrust, specific_impulse);
    intensity = Intensity(direction.coordinates() * speed);
  } else {
    CHECK(message.has_intensity()) << message;
    intensity = Intensity::ReadFromMessage(message.intensity());
  }

  Burn const burn{*intensity,
                  timing,
                  thrust,
                  specific_impulse,
                  RigidReferenceFrame<InertialFrame, Frame>::ReadFromMessage(
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
         frame()->FrenetFrame(t, to_frame_at_t(degrees_of_freedom))
             .template Forget<OrthogonalMap>();
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
typename Manœuvre<InertialFrame, Frame>::Timing const&
Manœuvre<InertialFrame, Frame>::full_timing() const {
  return burn_.timing;
}

template<typename InertialFrame, typename Frame>
Speed Manœuvre<InertialFrame, Frame>::ComputeЦиолковскийSpeed(
    Mass const& initial_mass,
    Time const& duration,
    Force const& thrust,
    SpecificImpulse const& specific_impulse) {
  Variation<Mass> const mass_flow = thrust / specific_impulse;
  Mass const final_mass = initial_mass - mass_flow * duration;
  return specific_impulse * std::log(initial_mass / final_mass);
}

}  // namespace internal
}  // namespace _manœuvre
}  // namespace ksp_plugin
}  // namespace principia
