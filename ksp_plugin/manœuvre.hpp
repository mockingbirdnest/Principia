#pragma once

#include "geometry/named_quantities.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Instant;
using geometry::Vector;
using physics::DiscreteTrajectory;
using quantities::Force;
using quantities::Mass;
using quantities::Speed;
using quantities::Time;
using quantities::Variation;

namespace ksp_plugin {

// This class represents a constant-thrust inertial burn.
template<typename Frame>
class Manœuvre {
 public:
  Manœuvre(Force thrust,
           Mass initial_mass,
           Speed effective_exhaust_velocity,
           Vector<double, Frame> direction);
  ~Manœuvre() = default;

  Force thrust() const;
  Mass initial_mass() const;
  // Specific impulse by mass, because specific impulse by weight is insane.
  // This is defined as the ratio of thrust to mass flow.
  // If the burn is done with a single engine (in a vacuum), this will be its
  // exhaust velocity.  For several engines, this is the total thrust divided
  // by the sum of the individual mass flows (where each mass flow is the
  // individual thrust divided by the exhaust velocity).
  Speed effective_exhaust_velocity() const;
  Vector<double, Frame> direction() const;

  // Equivalent characterizations of intensity.  Only one of the mutators may be
  // called, and only once.
  Time duration() const;
  void set_duration(Time const& duration);
  Speed Δv() const;
  void set_Δv(Speed const& Δv);

  // Equivalent characterizations of timing.  Only one of the mutators may be
  // called, and only once.
  Instant initial_time() const;
  void set_initial_time(Instant const& initial_time);
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
  Instant time_of_half_Δv() const;

  // Intensity and timing must have been set.
  typename DiscreteTrajectory<Frame>::IntrinsicAcceleration
  acceleration() const;

 private:
  Vector<double, Frame> const direction_;
  Speed const effective_exhaust_velocity_;
  std::unique_ptr<Time const> duration_;  // optional.
  std::unique_ptr<Instant const> initial_time_;  // optional.
  Force const thrust_;
  Mass const initial_mass_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/manœuvre_body.hpp"
