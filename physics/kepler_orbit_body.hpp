#pragma once

#include "physics/kepler_orbit.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {

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
  double const eccentricity = elements_at_epoch_.eccentricity;
  AngularFrequency const mean_motion =
      Sqrt(system_gravitational_parameter_ /
           Pow<3>(elements_at_epoch_.semimajor_axis)) * Radian;
  Angle const mean_anomaly = elements_at_epoch_.mean_nomaly +
                             mean_motion * (t - epoch_);
  auto const kepler_equation =
      [eccentricity, mean_anomaly](Angle const& eccentric_anomaly) -> Angle {
        return mean_anomaly -
                   (eccentric_anomaly - eccentricity * Sin(eccentric_anomaly));
      };
  Angle const eccentric_anomaly = FindRoot(kepler_equation);
  Angle const true_anomaly =
      ArcTan(Sqrt(1 - eccentricity) * Cos(eccentric_anomaly / 2),
             Sqrt(1 + eccentricity) * Sin(eccentric_anomaly / 2));
  
}

template<typename Argument, typename Function>
Argument FindRoot(Function f);

}  // namespace physics
}  // namespace principia
