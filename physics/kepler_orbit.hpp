#pragma once

#include "experimental/optional"

#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {


template<typename Frame>
struct KeplerianElements {
  double eccentricity;
  std::experimental::optional<Length> semimajor_axis;
  std::experimental::optional<AngularFrequency> mean_motion;
  Angle inclination;
  Angle longitude_of_ascending_node;
  Angle argument_of_periapsis;
  Angle mean_anomaly;
};

template<typename Frame>
class KeplerOrbit {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  // Exactly one of the |optional|s must be filled in  the given
  // |KeplerianElements|.
  KeplerOrbit(MassiveBody const& primary,
              Body const& secondary,
              Instant const& epoch,
              KeplerianElements<Frame> const& elements_at_epoch);

  // The |DegreesOfFreedom| of the secondary minus those of the primary.
  RelativeDegreesOfFreedom<Frame> StateVectors(Instant const& t) const;

  // All |optional|s are filled in the result.
  KeplerianElements<Frame> const& elements_at_epoch() const;

 private:
  GravitationalParameter gravitational_parameter_;
  KeplerianElements<Frame> elements_at_epoch_;
  Instant epoch_;
};

}  // namespace physics
}  // namespace principia

#include "physics/kepler_orbit_body.hpp"
