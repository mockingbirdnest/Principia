#pragma once

#include "physics/kepler_orbit.hpp"

#include "geometry/rotation.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {

using geometry::Bivector;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Time;

namespace physics {

template<typename Frame>
KeplerOrbit<Frame>::KeplerOrbit(
    TwoBodySystem const& system,
    Instant const& epoch,
    KeplerianElements<Frame> const& elements_at_epoch)
    : system_gravitational_parameter_(
          system.primary->gravitational_parameter() +
          system.secondary->is_massless()
              ? GravitationalParameter{}
              : CHECK_NOTNULL(
                    dynamic_cast<MassiveBody const*>(&*system.secondary))->
                        gravitational_parameter()),
      epoch_(epoch),
      elements_at_epoch_(element_at_epoch) {}

template<typename Frame>
KeplerOrbit<Frame>::StateVectors(Instant const& t) const {
  GravitationalParameter const μ = system_gravitational_parameter_;
  double const eccentricity = elements_at_epoch_.eccentricity;
  Length const a = elements_at_epoch_.semimajor_axis;
  Angle const i = elements_at_epoch_.inclination;
  Angle const Ω = elements_at_epoch_.longitude_of_ascending_node;
  Angle const ω = elements_at_epoch_.argument_of_periapsis;
  AngularFrequency const mean_motion = Sqrt(μ / Pow<3>(a)) * Radian;
  Angle const mean_anomaly = a + mean_motion * (t - epoch_);
  auto const kepler_equation =
      [eccentricity, mean_anomaly](Angle const& eccentric_anomaly) -> Angle {
        return mean_anomaly -
                   (eccentric_anomaly - eccentricity * Sin(eccentric_anomaly));
      };
  Angle const eccentric_anomaly = FindRoot(kepler_equation);
  // True anomaly.
  Angle const ν = ArcTan(Sqrt(1 - eccentricity) * Cos(eccentric_anomaly / 2),
                         Sqrt(1 + eccentricity) * Sin(eccentric_anomaly / 2));
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
      distance * from_orbit_plane(Vector<double, Frame>({Cos(ν), Sin(ν), 0}));
  Velocity<Frame> const v =
    Sqrt(μ * a) / distance *
        from_orbit_plane(
              Vector<double, Frame>(
                  {-Sin(eccentric_anomaly),
                   Sqrt(1 - Pow<2>(eccentricity)) * Cos(eccentric_anomaly),
                   0}));
  return {r, v};
}

// TODO(egg): move to numerics
template<typename Argument, typename Function>
Argument FindRoot(Function f, Argument first_guess) {

}

}  // namespace physics
}  // namespace principia
