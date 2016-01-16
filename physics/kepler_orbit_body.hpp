#pragma once

#include "physics/kepler_orbit.hpp"

#include "geometry/rotation.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {

using geometry::Bivector;
using numerics::Bisect;
using quantities::Cbrt;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Time;

namespace physics {

template<typename Frame>
KeplerOrbit<Frame>::KeplerOrbit(
    MassiveBody const& primary,
    Body const& secondary,
    KeplerianElements<Frame> const& elements_at_epoch,
    Instant const& epoch)
    : gravitational_parameter_(
          primary.gravitational_parameter() +
          (secondary.is_massless()
               ? GravitationalParameter{}
               : dynamic_cast<MassiveBody const&>(secondary).
                     gravitational_parameter())),
      elements_at_epoch_(elements_at_epoch),
      epoch_(epoch) {
  CHECK(static_cast<bool>(elements_at_epoch_.semimajor_axis) ^
        static_cast<bool>(elements_at_epoch_.mean_motion));
  GravitationalParameter const μ = gravitational_parameter_;
  if (elements_at_epoch_.semimajor_axis) {
    Length const& a = *elements_at_epoch_.semimajor_axis;
    elements_at_epoch_.mean_motion = Sqrt(μ / Pow<3>(a)) * Radian;
  } else {
    AngularFrequency const& n = *elements_at_epoch_.mean_motion;
    elements_at_epoch_.semimajor_axis = Cbrt(μ / Pow<2>(n / Radian));
  }
}

template<typename Frame>
RelativeDegreesOfFreedom<Frame>
KeplerOrbit<Frame>::StateVectors(Instant const& t) const {
  GravitationalParameter const& μ = gravitational_parameter_;
  double const& eccentricity = elements_at_epoch_.eccentricity;
  Length const& a = *elements_at_epoch_.semimajor_axis;
  Angle const& i = elements_at_epoch_.inclination;
  Angle const& Ω = elements_at_epoch_.longitude_of_ascending_node;
  Angle const& ω = elements_at_epoch_.argument_of_periapsis;
  Angle const mean_anomaly =
      elements_at_epoch_.mean_anomaly +
      *elements_at_epoch_.mean_motion * (t - epoch_);
  if (eccentricity < 1) {
    // Elliptic case.
    auto const kepler_equation =
        [eccentricity, mean_anomaly](Angle const& eccentric_anomaly) -> Angle {
          return mean_anomaly -
                     (eccentric_anomaly -
                      eccentricity * Sin(eccentric_anomaly) * Radian);
        };
    Angle const eccentric_anomaly =
        eccentricity == 0
            ? mean_anomaly
            : Bisect(kepler_equation,
                     mean_anomaly - eccentricity * Radian,
                     mean_anomaly + eccentricity * Radian);
    Angle const true_anomaly =
       2 * ArcTan(Sqrt(1 + eccentricity) * Sin(eccentric_anomaly / 2),
                  Sqrt(1 - eccentricity) * Cos(eccentric_anomaly / 2));
    Bivector<double, Frame> const x({1, 0, 0});
    Bivector<double, Frame> const y({0, 1, 0});
    Bivector<double, Frame> const z({0, 0, 1});
    // It would be nice to have a local frame, rather than make this a rotation
    // Frame -> Frame.
    // TODO(egg): Constructor for |Rotation| using Euler angles.
    Rotation<Frame, Frame> const from_orbit_plane =
        (Rotation<Frame, Frame>(Ω, z) *
         Rotation<Frame, Frame>(i, x) *
         Rotation<Frame, Frame>(ω, z));
    Length const distance = a * (1 - eccentricity * Cos(eccentric_anomaly));
    Displacement<Frame> const r =
        distance * from_orbit_plane(Vector<double, Frame>({Cos(true_anomaly),
                                                           Sin(true_anomaly),
                                                           0}));
    Velocity<Frame> const v =
        Sqrt(μ * a) / distance *
        from_orbit_plane(Vector<double, Frame>(
            {-Sin(eccentric_anomaly),
             Sqrt(1 - Pow<2>(eccentricity)) * Cos(eccentric_anomaly),
             0}));
    return {r, v};
  } else if (eccentricity == 1) {
    // Parabolic case.
    LOG(FATAL) << "not yet implemented";
    base::noreturn();
  } else {
    // Hyperbolic case.
    LOG(FATAL) << "not yet implemented";
    base::noreturn();
  }
}

template<typename Frame>
KeplerianElements<Frame> const& KeplerOrbit<Frame>::elements_at_epoch() const {
  return elements_at_epoch_;
}

}  // namespace physics
}  // namespace principia
