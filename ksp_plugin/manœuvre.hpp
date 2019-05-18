
#pragma once

#include <optional>

#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "physics/degrees_of_freedom.hpp"
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
using geometry::Velocity;
using physics::DegreesOfFreedom;
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

// This class represents a constant-thrust burn.  |InertialFrame| is an
// underlying inertial reference frame, |Frame| is the reference frame used to
// compute the Frenet frame.
template<typename InertialFrame, typename Frame>
class Manœuvre {
 public:
  // Characterization of intensity.  All members for exactly one of the groups
  // must be supplied.  The |direction| and |Δv| are given in the Frenet frame
  // of the trajectory at the beginning of the burn.
  struct Intensity final {
    // Group 1.
    std::optional<Vector<double, Frenet<Frame>>> direction;
    std::optional<Time> duration;
    // Group 2.
    std::optional<Velocity<Frenet<Frame>>> Δv;
  };

  // Characterization of timing.  All members for exactly one of the groups
  // must be supplied.
  struct Timing final {
    // Group 1.
    std::optional<Instant> initial_time;
    // Group 2.
    std::optional<Instant> time_of_half_Δv;
  };

  // Complete description of a burn.
  struct Burn final {
    Intensity intensity;
    Timing timing;
    Force thrust;
    // Specific impulse by mass, because specific impulse by weight is insane.
    // This is defined as the ratio of thrust to mass flow.
    // If the burn is done with a single engine (in a vacuum), this will be its
    // exhaust velocity.  For several engines, this is the total thrust divided
    // by the sum of the individual mass flows (where each mass flow is the
    // individual thrust divided by the exhaust velocity).
    SpecificImpulse specific_impulse;
    // Defines the Frenet frame.
    not_null<std::shared_ptr<DynamicFrame<InertialFrame, Frame> const>> frame;
    // If true, the direction of the burn remains fixed in a nonrotating frame.
    // Otherwise, the direction of the burn remains fixed in the Frenet frame of
    // the trajectory.
    bool is_inertially_fixed;
  };

  // A manœuvre is a burn applied to a vessel specified by its initial mass and
  // attached to a coasting trajectory.
  Manœuvre(Mass const& initial_mass,
           Burn const& burn);
  virtual ~Manœuvre() = default;

  Manœuvre(const Manœuvre&) = default;
  Manœuvre& operator=(const Manœuvre&) = default;

  Mass const& initial_mass() const;

  // Return the data passed at construction, with the same optionals set.
  Intensity const& intensity() const;
  Timing const& timing() const;
  Burn const& burn() const;

  // Individual intensity and timing fields.
  Vector<double, Frenet<Frame>> const& direction() const;
  Time const& duration() const;
  Velocity<Frenet<Frame>> const& Δv() const;
  Instant const& initial_time() const;
  Instant const& time_of_half_Δv() const;

  // Individual burn fields.
  Force const& thrust() const;
  SpecificImpulse const& specific_impulse() const;
  not_null<std::shared_ptr<DynamicFrame<InertialFrame, Frame> const>> frame()
      const;
  bool is_inertially_fixed() const;

  // Derived quantities.
  Variation<Mass> mass_flow() const;
  Mass final_mass() const;
  Time time_to_half_Δv() const;
  Instant final_time() const;

  // Returns true if ‖Δv‖ is NaN or infinite.
  bool IsSingular() const;

  // Returns true if and only if [initial_time, final_time] ⊆ ]begin, end[.
  bool FitsBetween(Instant const& begin, Instant const& end) const;

  // Sets the trajectory at the end of which the manœuvre takes place.  Must be
  // called before any of the functions below.  |trajectory| must have a point
  // at |initial_time()|.
  void set_coasting_trajectory(
      not_null<DiscreteTrajectory<InertialFrame> const*> trajectory);

  // This manœuvre must be inertially fixed.
  virtual Vector<double, InertialFrame> InertialDirection() const;

  // The result is valid until |*this| is destroyed.  This manœuvre must be
  // inertially fixed.
  typename Ephemeris<InertialFrame>::IntrinsicAcceleration
  InertialIntrinsicAcceleration() const;

  // Same as above for a manœuvre that is not inertially fixed.
  typename Ephemeris<InertialFrame>::GeneralizedIntrinsicAcceleration
  FrenetIntrinsicAcceleration() const;

  // Frenet frame at the beginning of the manœuvre.
  virtual OrthogonalMap<Frenet<Frame>, InertialFrame> FrenetFrame() const;

  // Note that |coasting_trajectory| is neither written nor read.
  void WriteToMessage(not_null<serialization::Manoeuvre*> message) const;
  static Manœuvre ReadFromMessage(
      serialization::Manoeuvre const& message,
      not_null<Ephemeris<InertialFrame>*> ephemeris);

 private:
  // Computes the Frenet frame at instant |t|, assuming that the motion has the
  // given |position| and |velocity|.
  OrthogonalMap<Frenet<Frame>, InertialFrame>
  ComputeFrenetFrame(
      Instant const& t,
      DegreesOfFreedom<InertialFrame> const& degrees_of_freedom) const;

  // Computes the acceleration at instant |t|, assuming that it happens in the
  // given |direction|.
  Vector<Acceleration, InertialFrame>
  ComputeIntrinsicAcceleration(
      Instant const& t,
      Vector<double, InertialFrame> const& direction) const;

  // Return structs where all the optionals are set.
  Intensity const& full_intensity() const;
  Timing const& full_timing() const;

  Mass initial_mass_;
  Burn construction_burn_;  // As given at construction.
  Burn burn_;  // All optionals filled.
  DiscreteTrajectory<InertialFrame> const* coasting_trajectory_ = nullptr;
};

}  // namespace internal_manœuvre

using internal_manœuvre::Manœuvre;

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/manœuvre_body.hpp"
