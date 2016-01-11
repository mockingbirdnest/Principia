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

  RelativeDegreesOfFreedom<Frame> StateVectors(Instant const& t) const;

 private:
  GravitationalParameter const system_gravitational_parameter_;
  KeplerianElements<Frame> const elements_at_epoch_;
  Instant const epoch_;
};

}  // namespace physics
}  // namespace principia

#include "physics/kepler_orbit_body.hpp"
