
#pragma once

#include <optional>

#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_manœuvre {

using base::not_null;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Vector;
using physics::DiscreteTrajectory;
using physics::DynamicFrame;
using physics::Ephemeris;
using physics::Frenet;
using quantities::Acceleration;
using quantities::Force;
using quantities::Mass;
using quantities::SpecificImpulse;
using quantities::Speed;
using quantities::Time;
using quantities::Variation;

// This class represents a constant-thrust inertial burn.  |InertialFrame| is
// an underlying inertial reference frame, |Frame| is the reference frame used
// to compute the Frenet frame.  |Frame| is defined by the parameter |frame|
// given to the constructor.  The |direction| is given in the Frenet frame of
// the trajectory at the beginning of the burn.
template<typename InertialFrame, typename Frame>
class Manœuvre {
 public:
  Manœuvre(Force const& thrust,
           Mass const& initial_mass,
           SpecificImpulse const& specific_impulse,
           Vector<double, Frenet<Frame>> const& direction,
           not_null<std::unique_ptr<DynamicFrame<InertialFrame, Frame> const>>
               frame,
           bool is_inertially_fixed);
  Manœuvre(Manœuvre&&) = default;
  Manœuvre& operator=(Manœuvre&&) = default;
  virtual ~Manœuvre() = default;

  Force const& thrust() const;
  Mass const& initial_mass() const;
  // Specific impulse by mass, because specific impulse by weight is insane.
  // This is defined as the ratio of thrust to mass flow.
  // If the burn is done with a single engine (in a vacuum), this will be its
  // exhaust velocity.  For several engines, this is the total thrust divided
  // by the sum of the individual mass flows (where each mass flow is the
  // individual thrust divided by the exhaust velocity).
  SpecificImpulse const& specific_impulse() const;
  Vector<double, Frenet<Frame>> const& direction() const;
  not_null<DynamicFrame<InertialFrame, Frame> const*> frame() const;

  // Equivalent characterizations of intensity.  Only one of the mutators may be
  // called, and only once.
  Time duration() const;
  void set_duration(Time const& duration);
  Speed Δv() const;
  void set_Δv(Speed const& Δv);

  // Intensity must have been set.  Returns true if Δv is NaN or infinite.
  bool IsSingular() const;

  // Equivalent characterizations of timing.  Only one of the mutators may be
  // called, and only once.
  Instant initial_time() const;
  void set_initial_time(Instant const& initial_time);
  // Intensity and timing must have been set.
  Instant time_of_half_Δv() const;
  // |Δv| or |duration| must have been set.
  void set_time_of_half_Δv(Instant const& time_of_half_Δv);

  // Derived quantities.
  Variation<Mass> mass_flow() const;
  // Intensity must have been set.
  Mass final_mass() const;
  // Intensity must have been set.
  Time time_to_half_Δv() const;
  // Intensity and timing must have been set.
  Instant final_time() const;

  // Intensity and timing must have been set.
  // Returns true if and only if [initial_time, final_time] ⊆ ]begin, end[.
  bool FitsBetween(Instant const& begin, Instant const& end) const;

  // Sets the trajectory at the end of which the mannœuvre takes place.  Must
  // be called before any of the functions below.  |trajectory| must have a
  // point at |initial_time()|.
  void set_coasting_trajectory(
      not_null<DiscreteTrajectory<InertialFrame> const*> trajectory);

  // If true, the direction of the burn remains fixed in a nonrotating frame.
  // Otherwise, the direction of the burn remains fixed in the Frenet frame of
  // the trajectory.
  bool is_inertially_fixed() const;

  // Intensity, timing and coasting trajectory must have been set.  This
  // manœuvre must be inertially fixed.
  virtual Vector<double, InertialFrame> InertialDirection() const;

  // Frenet frame at the beginning of the manœuvre.
  virtual OrthogonalMap<Frenet<Frame>, InertialFrame> FrenetFrame() const;

  // Intensity, timing and coasting trajectory must have been set.  The result
  // is valid until |*this| is destroyed.  This manœuvre must be inertially
  // fixed.
  typename Ephemeris<InertialFrame>::IntrinsicAcceleration
  IntrinsicAcceleration() const;

  // Intensity and timing must have been set.  |coasting_trajectory| is neither
  // written nor read.
  void WriteToMessage(not_null<serialization::Manoeuvre*> message) const;
  static Manœuvre ReadFromMessage(
      serialization::Manoeuvre const& message,
      not_null<Ephemeris<InertialFrame>*> ephemeris);

 private:
  // Computes the acceleration at |t|, assuming that it happens in the given
  // |direction|.
  Vector<Acceleration, InertialFrame>
  ComputeIntrinsicAcceleration(
      Instant const& t,
      Vector<double, InertialFrame> const& direction) const;

  Force const thrust_;
  Mass const initial_mass_;
  SpecificImpulse const specific_impulse_;
  Vector<double, Frenet<Frame>> const direction_;
  std::optional<Time> duration_;
  std::optional<Instant> initial_time_;
  not_null<std::unique_ptr<DynamicFrame<InertialFrame, Frame> const>> frame_;
  bool const is_inertially_fixed_;
  DiscreteTrajectory<InertialFrame> const* coasting_trajectory_ = nullptr;
};

}  // namespace internal_manœuvre

using internal_manœuvre::Manœuvre;

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/manœuvre_body.hpp"
