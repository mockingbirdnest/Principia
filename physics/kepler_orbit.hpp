#pragma once

#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {

template<typename Frame>
struct KeplerianElements {
  double eccentricity;
  Length semimajor_axis;
  Angle inclination;
  Angle longitude_of_ascending_node;
  Angle argument_of_periapsis;
  Angle mean_anomaly;
};

template<typename Frame>
class KeplerOrbit {
 public:
  KeplerOrbit(MassiveBody const& primary,
              Body const& secondary,
              Instant const& epoch,
              KeplerianElements<Frame> const& elements_at_epoch);

  RelativeDegreesOfFreedom<Frame> PrimocentricStateVectors(
      Instant const& t) const;
  RelativeDegreesOfFreedom<Frame> BarycentricStateVectors(
      Instant const& t) const;

 private:
  // The state vectors of a massless body orbiting a body with the given
  // |gravitation_parameter| in an orbit with the given |elements|.
  static RelativeDegreesOfFreedom<Frame> TestParticleStateVectors(
      KeplerianElements<Frame> const& elements,
      GravitationalParameter const& gravitational_parameter);
  GravitationalParameter const primary_gravitational_parameter_;
  GravitationalParameter const secondary_gravitational_parameter_;
  KeplerianElements<Frame> const elements_at_epoch_;
  Instant const epoch_;
};

}  // namespace physics
}  // namespace principia

#include "physics/kepler_orbit_body.hpp"
