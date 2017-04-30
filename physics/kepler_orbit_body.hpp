
#pragma once

#include "physics/kepler_orbit.hpp"

#include <string>

#include "geometry/rotation.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace physics {
namespace internal_kepler_orbit {

using geometry::AngleBetween;
using geometry::Bivector;
using geometry::Commutator;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::EulerAngles;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::OrientedAngleBetween;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using geometry::Wedge;
using numerics::Bisect;
using quantities::Abs;
using quantities::ArcCos;
using quantities::ArcSin;
using quantities::ArcTan;
using quantities::Cbrt;
using quantities::DebugString;
using quantities::Pow;
using quantities::Sin;
using quantities::SpecificAngularMomentum;
using quantities::SpecificEnergy;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Radian;

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
  GravitationalParameter const& μ = gravitational_parameter_;
  auto& eccentricity = elements_at_epoch_.eccentricity;
  auto& asymptotic_true_anomaly = elements_at_epoch_.asymptotic_true_anomaly;
  auto& turning_angle = elements_at_epoch_.turning_angle;
  auto& semimajor_axis = elements_at_epoch_.semimajor_axis;
  auto& specific_energy = elements_at_epoch_.specific_energy;
  auto& characteristic_energy = elements_at_epoch_.characteristic_energy;
  auto& mean_motion = elements_at_epoch_.mean_motion;
  auto& period = elements_at_epoch_.period;
  auto& hyperbolic_mean_motion = elements_at_epoch_.hyperbolic_mean_motion;
  auto& hyperbolic_excess_velocity =
      elements_at_epoch_.hyperbolic_excess_velocity;
  auto& semiminor_axis = elements_at_epoch_.semiminor_axis;
  auto& semilatus_rectum = elements_at_epoch_.semilatus_rectum;
  auto& specific_angular_momentum = elements_at_epoch_.specific_angular_momentum;
  auto& periapsis_distance = elements_at_epoch_.periapsis_distance;
  auto& apoapsis_distance = elements_at_epoch_.apoapsis_distance;
  int const eccentricity_specifications = eccentricity.has_value() +
                                          asymptotic_true_anomaly.has_value() +
                                          turning_angle.has_value();
  int const semimajor_axis_specifications =
      semimajor_axis.has_value() + specific_energy.has_value() +
      characteristic_energy.has_value() + mean_motion.has_value() +
      period.has_value() + hyperbolic_mean_motion.has_value() +
      hyperbolic_excess_velocity.has_value();
  int const semiminor_axis_specifications = semiminor_axis.has_value();
  int const semilatus_rectum_specifications =
      semilatus_rectum.has_value() + specific_angular_momentum.has_value();
  int const periapsis_distance_specifications = periapsis_distance.has_value();
  int const apoapsis_distance_specifications = apoapsis_distance.has_value();
  CHECK_LE(eccentricity_specifications, 1);
  CHECK_LE(semimajor_axis_specifications, 1);
  CHECK_LE(semiminor_axis_specifications, 1);
  CHECK_LE(semilatus_rectum_specifications, 1);
  CHECK_LE(periapsis_distance_specifications, 1);
  CHECK_LE(apoapsis_distance_specifications, 1);
  CHECK_EQ(eccentricity_specifications + semimajor_axis_specifications +
           semiminor_axis_specifications + semilatus_rectum_specifications +
           periapsis_distance_specifications + apoapsis_distance_specifications,
           2);
  // The conic shape and size is neither over- nor underspecified.
  if (eccentricity_specifications) {
    if (eccentricity) {
      double const& e = *eccentricity;
      turning_angle = 2 * ArcSin(1 / e);
      asymptotic_true_anomaly = ArcCos(1 / e);
    }
    if (turning_angle) {
      Angle const& δ = *turning_angle;
      eccentricity = 1 / Sin(δ / 2);
      asymptotic_true_anomaly = δ / 2 + π * Radian;
    }
    if (asymptotic_true_anomaly) {
      Angle const& θ_inf = *asymptotic_true_anomaly;
      eccentricity = -1 / Cos(θ_inf);
      turning_angle = 2 * (θ_inf - π * Radian);
    }
  }
  // TODO(egg): range checks.  What do we do with normalizable oddities
  // (negative n, T)? What about parabolae? (n = 0, a infinite, and the
  // equivalents provide an eccentricity specification...).
  if (semimajor_axis_specifications) {
    if (semimajor_axis) {
      Length const& a = *semimajor_axis;
      specific_energy = -μ / (2 * a);
      characteristic_energy = -μ / a;
      mean_motion = Sqrt(μ / Pow<3>(a)) * Radian;
      period = 2 * π * Sqrt(Pow<3>(a) / μ);
      hyperbolic_mean_motion = Sqrt(μ / Pow<3>(-a)) * Radian;
      hyperbolic_excess_velocity = Sqrt(-μ / a);
    }
    if (specific_energy) {
      SpecificEnergy const& ε = *specific_energy;
      semimajor_axis = -μ / (2 * ε);
      characteristic_energy = 2 * ε;
      mean_motion = 2 * Sqrt(-2 * Pow<3>(ε) / Pow<2>(μ)) * Radian;
      period = π * Sqrt(-Pow<2>(μ) / (2 * Pow<3>(ε)));
      hyperbolic_mean_motion = 2 * Sqrt(2 * Pow<3>(ε) / Pow<2>(μ)) * Radian;
      hyperbolic_excess_velocity = Sqrt(2 * ε);
    }
    if (characteristic_energy) {
      SpecificEnergy const& c3 = *characteristic_energy;
      semimajor_axis = -μ / c3;
      specific_energy = c3 / 2;
      mean_motion = Sqrt(-Pow<3>(c3) / Pow<2>(μ)) * Radian;
      period = 2 * π * Sqrt(-Pow<2>(μ) / Pow<3>(c3));
      hyperbolic_mean_motion = Sqrt(Pow<3>(c3) / Pow<2>(μ)) * Radian;
      hyperbolic_excess_velocity = Sqrt(c3);
    }
    if (mean_motion) {
      AngularFrequency const& n = *mean_motion;
      semimajor_axis = Cbrt(μ / Pow<2>(n / Radian));
      specific_energy = -Pow<2>(Cbrt(μ * n / Radian)) / 2;
      characteristic_energy =  -Pow<2>(Cbrt(μ * n / Radian));
      period = 2 * π * Radian / n;
      // The following two are NaN.
      hyperbolic_mean_motion = Sqrt(-Pow<2>(n));
      hyperbolic_excess_velocity = Sqrt(-Pow<2>(Cbrt(μ * n / Radian)));
    }
    if (period) {
      Time const& T = *period;
      semimajor_axis = Cbrt(μ * Pow<2>(T / (2 * π)));
      specific_energy = -Cbrt(Pow<2>(π * μ / T) / 2);
      characteristic_energy = -Cbrt(Pow<2>(2 * π * μ / T));
      mean_motion = 2 * π * Radian / T;
      // The following two are NaN.
      hyperbolic_mean_motion = 2 * π * Radian * Sqrt(-1 / Pow<2>(T));
      hyperbolic_excess_velocity = Cbrt(2 * π * Sqrt(-Pow<2>(μ / T)));
    }
    if (hyperbolic_mean_motion) {
      AngularFrequency const& n_over_i = *hyperbolic_mean_motion;
      semimajor_axis = -Cbrt(μ / Pow<2>(n_over_i / Radian));
      specific_energy = Pow<2>(Cbrt(μ * n_over_i / Radian)) / 2;
      characteristic_energy =  Pow<2>(Cbrt(μ * n_over_i / Radian));
      // The following two are NaN.
      period = 2 * π * Radian / Sqrt(-Pow<2>(n_over_i));
      mean_motion = Sqrt(-Pow<2>(n_over_i));
      hyperbolic_excess_velocity = Sqrt(Pow<2>(Cbrt(μ * n_over_i / Radian)));
    }
    if (hyperbolic_excess_velocity) {
      Speed const& v_inf = *hyperbolic_excess_velocity;
      semimajor_axis = -μ / Pow<2>(v_inf);
      specific_energy = Pow<2>(v_inf) / 2;
      characteristic_energy =  Pow<2>(v_inf);
      // The following two are NaN.
      period = 2 * π * Sqrt(-Pow<2>(μ) / Pow<6>(v_inf));
      mean_motion = Sqrt(-Pow<6>(v_inf) / Pow<2>(μ)) * Radian;
      hyperbolic_mean_motion = Sqrt(Pow<6>(v_inf) / Pow<2>(μ)) * Radian;
    }
  }
  if (semilatus_rectum_specifications) {
    if (semilatus_rectum) {
      specific_angular_momentum = Sqrt(μ * *semilatus_rectum) / Radian;
    }
    if (specific_angular_momentum) {
      SpecificAngularMomentum const& h = *specific_angular_momentum;
      semilatus_rectum = Pow<2>(h * Radian) / μ;
    }
  }

  if (eccentricity_specifications && semimajor_axis_specifications) {
    double const& e = *eccentricity;
    Length const& a = *semimajor_axis;
    semiminor_axis = a * Sqrt(Abs(1 - Pow<2>(e)));
    semilatus_rectum = a * (1 - Pow<2>(e));
    periapsis_distance = a * (1 - e);
    apoapsis_distance = a * e;
  }
  if (eccentricity_specifications && semiminor_axis_specifications) {
    double const& e = *eccentricity;
    Length const& b = *semiminor_axis;
    semimajor_axis = b / Sqrt(Abs(1 - Pow<2>(e)));
    semilatus_rectum = Abs(b) * Sqrt(Abs(1 - Pow<2>(e)));
    periapsis_distance = *semimajor_axis * (1 - e);
    periapsis_distance = *semimajor_axis * e;
  }
  if (eccentricity_specifications && semilatus_rectum_specifications) {
    double const& e = *eccentricity;
    Length const& p = *semilatus_rectum;
    semimajor_axis = p / (1 - Pow<2>(e));
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

  Bivector<SpecificAngularMomentum, Frame> const h = Wedge(r / Radian, v);
  // The eccentricity vector has magnitude equal to the eccentricity, and points
  // towards the periapsis.  This is a vector (the direction of the periapsis
  // does not depend on the coordinate system).
  Vector<double, Frame> const eccentricity_vector =
      v * h / μ * Radian - Normalize(r);
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
    struct OrbitPlane;
    Rotation<OrbitPlane, Frame> const from_orbit_plane(
        Ω, i, ω,
        EulerAngles::ZXZ,
        DefinesFrame<OrbitPlane>{});
    Length const distance = a * (1 - eccentricity * Cos(eccentric_anomaly));
    Displacement<Frame> const r =
        distance * from_orbit_plane(Vector<double, OrbitPlane>(
                       {Cos(true_anomaly), Sin(true_anomaly), 0}));
    Velocity<Frame> const v =
        Sqrt(μ * a) / distance *
        from_orbit_plane(Vector<double, OrbitPlane>(
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

}  // namespace internal_kepler_orbit
}  // namespace physics
}  // namespace principia
