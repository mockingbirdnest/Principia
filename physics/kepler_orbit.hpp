
#pragma once

#include <optional>
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
  std::optional<double> eccentricity;
  // The following two elements are NaN for elliptic orbits.
  std::optional<Angle> asymptotic_true_anomaly;
  std::optional<Angle> turning_angle;
  // 2. semimajor axis.
  std::optional<Length> semimajor_axis;
  std::optional<SpecificEnergy> specific_energy;
  std::optional<SpecificEnergy> characteristic_energy;
  // The following two elements are NaN for hyperbolic orbits.
  std::optional<AngularFrequency> mean_motion;
  std::optional<Time> period;
  // The following two elements are NaN for elliptic orbits.
  std::optional<AngularFrequency> hyperbolic_mean_motion;
  std::optional<Speed> hyperbolic_excess_velocity;
  // 3. semiminor axis.  The |semiminor_axis| is NaN for hyperbolic orbits, the
  // |impact_parameter| is NaN for elliptic orbits.
  std::optional<Length> semiminor_axis;
  std::optional<Length> impact_parameter;
  // 4. semilatus rectum.
  std::optional<Length> semilatus_rectum;
  std::optional<SpecificAngularMomentum>
      specific_angular_momentum;
  // 5. periapsis distance.
  std::optional<Length> periapsis_distance;
  // 6. apoapsis distance.
  std::optional<Length> apoapsis_distance;

  // II. These elements determine the orientation of the conic.  Three are
  // needed.
  Angle inclination;
  Angle longitude_of_ascending_node;
  std::optional<Angle> argument_of_periapsis;
  std::optional<Angle> longitude_of_periapsis;

  // III. These elements determine a point on the conic.  One is needed.
  std::optional<Angle> true_anomaly;
#if NOT_YET_IMPLEMENTED
  std::optional<Angle> true_longitude;
  std::optional<Time> time_since_periapsis;
#endif
  // The mean anomaly and mean longitude are NaN for hyperbolic orbits.
  std::optional<Angle> mean_anomaly;
#if NOT_YET_IMPLEMENTED
  std::optional<Angle> mean_longitude;
#endif
  // The hyperbolic mean anomaly is NaN for elliptic orbits.
  std::optional<Angle> hyperbolic_mean_anomaly;

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
  // |elements| must be minimally specified.  Fills all |optional|s in
  // |elements|.
  static void CompleteElements(KeplerianElements<Frame>& elements,
                               GravitationalParameter const& μ);
  // For each category in section I of |elements|, either one, none, or all of
  // the |optional|s must be filled.  If one is filled, fills the others in that
  // category.
  static void CompleteConicParametersByCategory(
      KeplerianElements<Frame>& elements,
      GravitationalParameter const& μ);
  // Section I of |elements| must be minimally specified.  Fills it.
  static void CompleteConicParameters(KeplerianElements<Frame>& elements,
                                      GravitationalParameter const& μ);
  // Section II of |elements| must be minimally specified.  Fills it.
  static void CompleteOrientationParameters(KeplerianElements<Frame>& elements);
  // Sections I and II of |elements| must be filled; section III must be
  // minimally specified.  Fills section III.
  static void CompleteAnomalies(KeplerianElements<Frame>& elements);

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
