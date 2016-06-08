
#pragma once

#include "physics/kepler_orbit.hpp"

#include <string>

#include "geometry/rotation.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {

using geometry::AngleBetween;
using geometry::Bivector;
using geometry::Commutator;
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::OrientedAngleBetween;
using geometry::Wedge;
using numerics::Bisect;
using quantities::ArcCos;
using quantities::Cbrt;
using quantities::DebugString;
using quantities::Pow;
using quantities::SpecificAngularMomentum;
using quantities::SpecificEnergy;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Time;

namespace physics {

template<typename Frame>
void KeplerianElements<Frame>::WriteToMessage(
    not_null<serialization::KeplerianElements*> const message) const {
  Frame::WriteToMessage(message->mutable_frame());
  message->set_eccentricity(eccentricity);
  if (semimajor_axis) {
    semimajor_axis->WriteToMessage(message->mutable_semimajor_axis());
  }
  if (mean_motion) {
    mean_motion->WriteToMessage(message->mutable_mean_motion());
  }
  inclination.WriteToMessage(message->mutable_inclination());
  longitude_of_ascending_node.WriteToMessage(
      message->mutable_longitude_of_ascending_node());
  argument_of_periapsis.WriteToMessage(
      message->mutable_argument_of_periapsis());
  mean_anomaly.WriteToMessage(message->mutable_mean_anomaly());
}

template<typename Frame>
std::string DebugString(KeplerianElements<Frame> const& elements) {
  std::string result = "{";
  auto const append = [&result](
      std::string const& symbol,
      auto const& value,
      bool end = false) {
    result += symbol + " = " + quantities::DebugString(value) +
              (end ? "}" : ", ");
  };
  append("e", elements.eccentricity);
  if (elements.semimajor_axis) {
    append("a", *elements.semimajor_axis);
  }
  if (elements.mean_motion) {
    append("n", *elements.mean_motion);
  }
  append("i", elements.inclination);
  append(u8"Ω", elements.longitude_of_ascending_node);
  append(u8"ω", elements.argument_of_periapsis);
  append("M", elements.mean_anomaly, /*end=*/true);
  return result;
}

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         KeplerianElements<Frame> const& elements) {
  return out << DebugString(elements);
}

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
KeplerOrbit<Frame>::KeplerOrbit(
    MassiveBody const& primary,
    Body const& secondary,
    RelativeDegreesOfFreedom<Frame> const& state_vectors,
    Instant const& epoch)
    : gravitational_parameter_(
          primary.gravitational_parameter() +
          (secondary.is_massless() ? GravitationalParameter{}
                                   : dynamic_cast<MassiveBody const&>(secondary)
                                         .gravitational_parameter())),
      epoch_(epoch) {
  GravitationalParameter const& μ = gravitational_parameter_;
  Displacement<Frame> const& r = state_vectors.displacement();
  Velocity<Frame> const& v = state_vectors.velocity();
  Vector<double, Frame> const x({1, 0, 0});
  Vector<double, Frame> const z({0, 0, 1});
  Bivector<double, Frame> const x_wedge_y({0, 0, 1});

  Bivector<SpecificAngularMomentum, Frame> const h = Wedge(r, v) * Radian;
  // The eccentricity vector has magnitude equal to the eccentricity, and points
  // towards the periapsis.  This is a vector (the direction of the periapsis
  // does not depend on the coordinate system).
  Vector<double, Frame> const eccentricity_vector =
      v * h / μ / Radian - Normalize(r);
  auto const& periapsis = eccentricity_vector;
  Vector<SpecificAngularMomentum, Frame> const ascending_node = z * h;

  // Maps [-π, π] to [0, 2π].
  auto const positive_angle = [](Angle const& α) -> Angle {
    return α > 0 * Radian ? α : α + 2 * π * Radian;
  };

  // Inclination (above the xy plane).
  Angle const i = AngleBetween(x_wedge_y, h);
  // Argument of periapsis.
  Angle const ω = positive_angle(
      OrientedAngleBetween(ascending_node, periapsis, x_wedge_y));
  // Longitude of ascending node.
  // This is equivalent to |OrientedAngleBetween(x, ascending_node, x_wedge_y)|
  // since |ascending_node| lies in the xy plane.
  Angle const Ω = positive_angle(
      ArcTan(ascending_node.coordinates().y, ascending_node.coordinates().x));
  double const eccentricity = eccentricity_vector.Norm();
  Angle const true_anomaly =
      positive_angle(OrientedAngleBetween(periapsis, r, x_wedge_y));
  Angle const eccentric_anomaly =
      ArcTan(Sqrt(1 - Pow<2>(eccentricity)) * Sin(true_anomaly),
             eccentricity + Cos(true_anomaly));
  Angle const mean_anomaly = positive_angle(
      eccentric_anomaly - eccentricity * Sin(eccentric_anomaly) * Radian);

  SpecificEnergy const ε = InnerProduct(v, v) / 2 - μ / r.Norm();
  // Semimajor axis.
  Length const a = -μ / (2 * ε);
  // Mean motion.
  AngularFrequency const n = Sqrt(μ / Pow<3>(a)) * Radian;

  elements_at_epoch_.eccentricity                = eccentricity;
  elements_at_epoch_.semimajor_axis              = a;
  elements_at_epoch_.mean_motion                 = n;
  elements_at_epoch_.inclination                 = i;
  elements_at_epoch_.longitude_of_ascending_node = Ω;
  elements_at_epoch_.argument_of_periapsis       = ω;
  elements_at_epoch_.mean_anomaly                = mean_anomaly;
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
