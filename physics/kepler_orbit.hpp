
#pragma once

#include <experimental/optional>
#include <ostream>
#include <string>

#include "physics/body.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {
namespace internal_kepler_orbit {

using base::not_null;
using geometry::Instant;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::SpecificAngularMomentum;
using quantities::SpecificEnergy;
using quantities::Speed;
using quantities::Time;

template<typename Frame>
struct KeplerianElements final {
  // I. These elements determine the shape and size of the conic.  Two are
  // needed, from two different numbered categories below.
  // 1. eccentricity.
  std::experimental::optional<double> eccentricity;
  // The following two elements are NaN for elliptic orbits.
  std::experimental::optional<Angle> asymptotic_true_anomaly;
  std::experimental::optional<Angle> turning_angle;
  // 2. semimajor axis.
  std::experimental::optional<Length> semimajor_axis;
  std::experimental::optional<SpecificEnergy> specific_energy;
  std::experimental::optional<SpecificEnergy> characteristic_energy;
  // The following two elements are NaN for hyperbolic orbits.
  std::experimental::optional<AngularFrequency> mean_motion;
  std::experimental::optional<Time> period;
  // The following two elements are NaN for elliptic orbits.
  std::experimental::optional<AngularFrequency> hyperbolic_mean_motion;
  std::experimental::optional<Speed> hyperbolic_excess_velocity;
  // 3. semiminor axis.  The |semiminor_axis| is NaN for hyperbolic orbits, the
  // |impact_parameter| is NaN for elliptic orbits.
  std::experimental::optional<Length> semiminor_axis;
  std::experimental::optional<Length> impact_parameter;
  // 4. semilatus rectum.
  std::experimental::optional<Length> semilatus_rectum;
  std::experimental::optional<SpecificAngularMomentum>
      specific_angular_momentum;
  // 5. periapsis distance.
  std::experimental::optional<Length> periapsis_distance;
  // 6. apoapsis distance.
  std::experimental::optional<Length> apoapsis_distance;

  // II. These elements determine the orientation of the conic.  Three are
  // needed.
  Angle inclination;
  Angle longitude_of_ascending_node;
  std::experimental::optional<Angle> argument_of_periapsis;
  std::experimental::optional<Angle> longitude_of_periapsis;

  // III. These elements determine a point on the conic.  One is needed.
  std::experimental::optional<Angle> true_anomaly;
  std::experimental::optional<Angle> true_longitude;
  std::experimental::optional<Time> time_since_periapsis;
  // The mean anomaly and mean longitude are NaN for hyperbolic orbits.
  std::experimental::optional<Angle> mean_anomaly;
  std::experimental::optional<Angle> mean_longitude;
  // The hyperbolic mean anomaly is NaN for elliptic orbits:
  // |hyperbolic_mean_anomaly = i * mean_anomaly|.
  std::experimental::optional<Angle> hyperbolic_mean_anomaly;

  void WriteToMessage(
      not_null<serialization::KeplerianElements*> message) const;
};

template<typename Frame>
std::string DebugString(KeplerianElements<Frame> const& elements);

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         KeplerianElements<Frame> const& elements);

template<typename Frame>
class KeplerOrbit final {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  // Exactly one of the |optional|s must be filled in the given
  // |KeplerianElements|.
  KeplerOrbit(MassiveBody const& primary,
              Body const& secondary,
              KeplerianElements<Frame> const& elements_at_epoch,
              Instant const& epoch);
  KeplerOrbit(MassiveBody const& primary,
              Body const& secondary,
              RelativeDegreesOfFreedom<Frame> const& state_vectors,
              Instant const& epoch);

  // The |DegreesOfFreedom| of the secondary minus those of the primary.
  RelativeDegreesOfFreedom<Frame> StateVectors(Instant const& t) const;

  // All |optional|s are filled in the result.
  KeplerianElements<Frame> const& elements_at_epoch() const;

 private:
  GravitationalParameter const gravitational_parameter_;
  KeplerianElements<Frame> elements_at_epoch_;
  Instant const epoch_;
};

}  // namespace internal_kepler_orbit

using internal_kepler_orbit::KeplerianElements;
using internal_kepler_orbit::KeplerOrbit;

}  // namespace physics
}  // namespace principia

#include "physics/kepler_orbit_body.hpp"
