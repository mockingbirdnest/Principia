#pragma once

#include "geometry/named_quantities.hpp"
#include "physics/trajectory.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using geometry::Instant;
using geometry::Vector;
using physics::Trajectory;
using quantities::Force;
using quantities::Mass;
using quantities::Speed;
using quantities::Time;
using quantities::Variation;

namespace ksp_plugin {

// This class represents a constant-thrust inertial burn.
template<typename Frame>
class Manœuvre {
  // Equivalent characterizations of timing.
  Instant start_time() const;
  Instant half_Δv() const;
  Instant end_time() const;

  // Equivalent characterizations of intensity.
  Time duration() const;
  Speed Δv() const;

  Vector<double, Frame> direction() const;

  // Specific impulse by mass, because specific impulse by weight is insane.
  Speed exhaust_velocity() const;

  Force thrust() const;

  Mass initial_mass() const;

  // Derived quantities.
  Variation<Mass> mass_flow() const;
  Mass final_mass() const;
  double mass_ratio() const;

  typename Trajectory<Frame>::IntrinsicAcceleration acceleration() const;
 private:
  Vector<double, Frame> direction_;
  Instant start_time_;
  Instant end_time_;
  Speed exhaust_velocity_;
  Force thrust_;
  Mass initial_mass_;
};

}  // namespace ksp_plugin
}  // namespace principia
