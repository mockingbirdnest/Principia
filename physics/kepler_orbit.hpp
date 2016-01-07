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

struct TwoBodySystem {
  not_null<MassiveBody const*> const primary;
  not_null<Body const*> const secondary;
};

template<typename Frame>
class KeplerOrbit {
 public:
  KeplerOrbit(TwoBodySystem const& system,
              Instant const& epoch,
              KeplerianElements<Frame> const& elements_at_epoch);
  StateVectors(Instant const& t) const;
 private:
  GravitationalParameter const sysetm_gravitational_parameter_;
  KeplerianElements const elements_at_epoch_;
  Instant const epoch_;
};

}  // namespace physics
}  // namespace principia

#include "physics/kepler_orbit_body.hpp"
