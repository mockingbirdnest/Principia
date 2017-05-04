
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
using quantities::ArcCosh;
using quantities::ArcSin;
using quantities::ArcTan;
using quantities::Cbrt;
using quantities::DebugString;
using quantities::NaN;
using quantities::Pow;
using quantities::Sin;
using quantities::Sinh;
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
  bool first_entry = true;
  auto const append = [&result, &first_entry](std::string const& symbol,
                                              auto const& value) {
    result += (first_entry ? "" : ", ") + symbol + " = " +
              quantities::DebugString(value);
    first_entry = false;
  };
  auto const append_optional = [&append](std::string const& symbol,
                                         auto const& value) {
    if (value) {
      append(symbol, *value);
    }
  };
  append_optional("e", elements.eccentricity);
  append_optional(u8"ν∞", elements.asymptotic_true_anomaly);
  append_optional(u8"δ", elements.turning_angle);

  append_optional("a", elements.semimajor_axis);
  append_optional(u8"ε", elements.specific_energy);
  append_optional(u8"C₃", elements.characteristic_energy);
  append_optional("n", elements.mean_motion);
  append_optional("T", elements.period);
  append_optional("n/i", elements.hyperbolic_mean_motion);
  append_optional(u8"v∞", elements.hyperbolic_excess_velocity);

  append_optional("b", elements.semiminor_axis);
  append_optional("b/i", elements.impact_parameter);

  append_optional("p", elements.semilatus_rectum);
  append_optional("h", elements.specific_angular_momentum);

  append_optional("r_pe", elements.periapsis_distance);

  append_optional("r_ap", elements.apoapsis_distance);

  append("i", elements.inclination);
  append(u8"Ω", elements.longitude_of_ascending_node);
  append_optional(u8"ω", elements.argument_of_periapsis);
  append_optional(u8"ϖ", elements.longitude_of_periapsis);

  append_optional(u8"ν", elements.true_anomaly);
  append_optional("M", elements.mean_anomaly);
  append_optional("M/i", elements.hyperbolic_mean_anomaly);
  result += "}";
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
  CompleteElements(elements_at_epoch_, gravitational_parameter_);
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
  Angle const true_anomaly =
      positive_angle(OrientedAngleBetween(periapsis, r, x_wedge_y));

  SpecificEnergy const ε = InnerProduct(v, v) / 2 - μ / r.Norm();
  double e = eccentricity_vector.Norm();

  // We have h, e, and ε.  There are three ways of computing each of b, r_pe,
  // and r_ap from that (from any two of the elements we have), but only one is
  // well-conditioned for all eccentricities.
  Length const b = Sqrt(-InnerProduct(h, h) / (2 * ε));
  Length const impact_parameter = Sqrt(InnerProduct(h, h) / (2 * ε));
  Length const r_pe = InnerProduct(h, h) / ((1 + e) * μ);
  Length const r_ap = - μ * (1 + e) / (2 * ε);

  elements_at_epoch_.eccentricity                = e;
  elements_at_epoch_.specific_angular_momentum   = h.Norm();
  elements_at_epoch_.semiminor_axis              = b;
  elements_at_epoch_.impact_parameter            = impact_parameter;
  elements_at_epoch_.periapsis_distance          = r_pe;
  elements_at_epoch_.apoapsis_distance           = r_ap;
  elements_at_epoch_.inclination                 = i;
  elements_at_epoch_.longitude_of_ascending_node = Ω;
  elements_at_epoch_.argument_of_periapsis       = ω;
  elements_at_epoch_.true_anomaly                = true_anomaly;
  CompleteElements(elements_at_epoch_, μ);
}

template<typename Frame>
RelativeDegreesOfFreedom<Frame>
KeplerOrbit<Frame>::StateVectors(Instant const& t) const {
  GravitationalParameter const& μ = gravitational_parameter_;
  double const& e = *elements_at_epoch_.eccentricity;
  Angle const& i = elements_at_epoch_.inclination;
  Angle const& Ω = elements_at_epoch_.longitude_of_ascending_node;
  Angle const& ω = *elements_at_epoch_.argument_of_periapsis;
  Length const& p = *elements_at_epoch_.semilatus_rectum;
  SpecificEnergy const& ε = *elements_at_epoch_.specific_energy;
  KeplerianElements<Frame> elements = elements_at_epoch_;
  elements.true_anomaly.reset();
  elements.mean_anomaly.reset();
  elements.hyperbolic_mean_anomaly.reset();
  if (e < 1) {
    elements.mean_anomaly = *elements_at_epoch_.mean_anomaly +
                            *elements_at_epoch_.mean_motion * (t - epoch_);
  } else if (e == 1) {
    // Parabolic case.
    LOG(FATAL) << "not yet implemented";
  } else {
    elements.hyperbolic_mean_anomaly =
        *elements_at_epoch_.hyperbolic_mean_anomaly +
        *elements_at_epoch_.hyperbolic_mean_motion * (t - epoch_);
  }
  CompleteAnomalies(elements);
  Angle const& ν = *elements.true_anomaly;
  struct OrbitPlane;
  Rotation<OrbitPlane, Frame> const from_orbit_plane(
      Ω, i, ω,
      EulerAngles::ZXZ,
      DefinesFrame<OrbitPlane>{});
  Length const r = p / (1 + e * Cos(ν));
  Displacement<Frame> const displacement =
      r * from_orbit_plane(Vector<double, OrbitPlane>({Cos(ν), Sin(ν), 0}));
  // Flight path angle.
  Angle const φ = ArcTan(e * Sin(ν), 1 + e * Cos(ν));
  // The norm comes from the vis-viva equation.
  Velocity<Frame> const velocity =
      Sqrt(2 * (ε + μ / r)) *
      from_orbit_plane(Vector<double, OrbitPlane>(
          {-Sin(ν - φ), Cos(ν - φ), 0}));
  return {displacement, velocity};
}

template<typename Frame>
KeplerianElements<Frame> const& KeplerOrbit<Frame>::elements_at_epoch() const {
  return elements_at_epoch_;
}

template<typename Frame>
void KeplerOrbit<Frame>::CompleteElements(KeplerianElements<Frame>& elements,
                                          GravitationalParameter const& μ) {
  auto& eccentricity = elements.eccentricity;
  auto& asymptotic_true_anomaly = elements.asymptotic_true_anomaly;
  auto& turning_angle = elements.turning_angle;
  auto& semimajor_axis = elements.semimajor_axis;
  auto& specific_energy = elements.specific_energy;
  auto& characteristic_energy = elements.characteristic_energy;
  auto& mean_motion = elements.mean_motion;
  auto& period = elements.period;
  auto& hyperbolic_mean_motion = elements.hyperbolic_mean_motion;
  auto& hyperbolic_excess_velocity = elements.hyperbolic_excess_velocity;
  auto& semiminor_axis = elements.semiminor_axis;
  auto& impact_parameter = elements.impact_parameter;
  auto& semilatus_rectum = elements.semilatus_rectum;
  auto& specific_angular_momentum = elements.specific_angular_momentum;
  auto& periapsis_distance = elements.periapsis_distance;
  auto& apoapsis_distance = elements.apoapsis_distance;
  int const eccentricity_specifications = eccentricity.has_value() +
                                          asymptotic_true_anomaly.has_value() +
                                          turning_angle.has_value();
  int const semimajor_axis_specifications =
      semimajor_axis.has_value() + specific_energy.has_value() +
      characteristic_energy.has_value() + mean_motion.has_value() +
      period.has_value() + hyperbolic_mean_motion.has_value() +
      hyperbolic_excess_velocity.has_value();
  int const semiminor_axis_specifications =
      semiminor_axis.has_value() + impact_parameter.has_value();
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
  bool const eccentricity_specified = eccentricity_specifications;
  bool const semimajor_axis_specified = semimajor_axis_specifications;
  bool const semiminor_axis_specified = semiminor_axis_specifications;
  bool const semilatus_rectum_specified = semilatus_rectum_specifications;
  bool const periapsis_distance_specified = periapsis_distance_specifications;
  bool const apoapsis_distance_specified = apoapsis_distance_specifications;
  CHECK_EQ(eccentricity_specified + semimajor_axis_specified +
           semiminor_axis_specified + semilatus_rectum_specified +
           periapsis_distance_specified + apoapsis_distance_specified,
           2);
  // In the first pass, we fill all the parameters within the categories that
  // are specified.  We then specify one parameter each remaining category, and
  // in a second pass, we fill those categories.
  // NOTE(egg): implicit capture because we want all of the renamings and
  // booleans above.
  auto const complete_conic_parameters = [&](bool const first_pass) {
    // The conic shape and size is neither over- nor underspecified.
    if (eccentricity_specified == first_pass) {
      if (eccentricity) {
        double const& e = *eccentricity;
        turning_angle = 2 * ArcSin(1 / e);
        asymptotic_true_anomaly = ArcCos(1 / e);
      } else if (turning_angle) {
        Angle const& δ = *turning_angle;
        eccentricity = 1 / Sin(δ / 2);
        asymptotic_true_anomaly = δ / 2 + π * Radian;
      } else if (asymptotic_true_anomaly) {
        Angle const& ν_inf = *asymptotic_true_anomaly;
        eccentricity = -1 / Cos(ν_inf);
        turning_angle = 2 * (ν_inf - π * Radian);
      }
    }
    // TODO(egg): range checks.  What do we do with normalizable oddities
    // (negative n, T)? What about parabolae? (n = 0, a infinite, and the
    // equivalents provide an eccentricity specification...).
    if (semimajor_axis_specified == first_pass) {
      if (semimajor_axis) {
        Length const& a = *semimajor_axis;
        specific_energy = -μ / (2 * a);
        characteristic_energy = -μ / a;
        mean_motion = Sqrt(μ / Pow<3>(a)) * Radian;
        period = 2 * π * Sqrt(Pow<3>(a) / μ);
        hyperbolic_mean_motion = Sqrt(μ / Pow<3>(-a)) * Radian;
        hyperbolic_excess_velocity = Sqrt(-μ / a);
      } else if (specific_energy) {
        SpecificEnergy const& ε = *specific_energy;
        semimajor_axis = -μ / (2 * ε);
        characteristic_energy = 2 * ε;
        mean_motion = 2 * Sqrt(-2 * Pow<3>(ε) / Pow<2>(μ)) * Radian;
        period = π * Sqrt(-Pow<2>(μ) / (2 * Pow<3>(ε)));
        hyperbolic_mean_motion = 2 * Sqrt(2 * Pow<3>(ε) / Pow<2>(μ)) * Radian;
        hyperbolic_excess_velocity = Sqrt(2 * ε);
      } else if (characteristic_energy) {
        SpecificEnergy const& c3 = *characteristic_energy;
        semimajor_axis = -μ / c3;
        specific_energy = c3 / 2;
        mean_motion = Sqrt(-Pow<3>(c3) / Pow<2>(μ)) * Radian;
        period = 2 * π * Sqrt(-Pow<2>(μ) / Pow<3>(c3));
        hyperbolic_mean_motion = Sqrt(Pow<3>(c3) / Pow<2>(μ)) * Radian;
        hyperbolic_excess_velocity = Sqrt(c3);
      } else if (mean_motion) {
        AngularFrequency const& n = *mean_motion;
        semimajor_axis = Cbrt(μ / Pow<2>(n / Radian));
        specific_energy = -Pow<2>(Cbrt(μ * n / Radian)) / 2;
        characteristic_energy =  -Pow<2>(Cbrt(μ * n / Radian));
        period = 2 * π * Radian / n;
        // The following two are NaN.
        hyperbolic_mean_motion = Sqrt(-Pow<2>(n));
        hyperbolic_excess_velocity = Sqrt(-Pow<2>(Cbrt(μ * n / Radian)));
      } else if (period) {
        Time const& T = *period;
        semimajor_axis = Cbrt(μ * Pow<2>(T / (2 * π)));
        specific_energy = -Cbrt(Pow<2>(π * μ / T) / 2);
        characteristic_energy = -Cbrt(Pow<2>(2 * π * μ / T));
        mean_motion = 2 * π * Radian / T;
        // The following two are NaN.
        hyperbolic_mean_motion = 2 * π * Radian * Sqrt(-1 / Pow<2>(T));
        hyperbolic_excess_velocity = Cbrt(2 * π * Sqrt(-Pow<2>(μ / T)));
      } else if (hyperbolic_mean_motion) {
        AngularFrequency const& n_over_i = *hyperbolic_mean_motion;
        semimajor_axis = -Cbrt(μ / Pow<2>(n_over_i / Radian));
        specific_energy = Pow<2>(Cbrt(μ * n_over_i / Radian)) / 2;
        characteristic_energy =  Pow<2>(Cbrt(μ * n_over_i / Radian));
        // The following two are NaN.
        period = 2 * π * Radian / Sqrt(-Pow<2>(n_over_i));
        mean_motion = Sqrt(-Pow<2>(n_over_i));
        hyperbolic_excess_velocity = Sqrt(Pow<2>(Cbrt(μ * n_over_i / Radian)));
      } else if (hyperbolic_excess_velocity) {
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
    // Because they differ by a factor of i, we fill both of those between the
    // first and the second pass, thus this is only needed in the first pass.
    if(semiminor_axis_specified && first_pass) {
      if (semiminor_axis) {
        impact_parameter = Sqrt(-Pow<2>(*semiminor_axis));
      } else if (impact_parameter) {
        semiminor_axis = Sqrt(-Pow<2>(*impact_parameter));
      }
    }
    if (semilatus_rectum_specified == first_pass) {
      if (semilatus_rectum) {
        specific_angular_momentum = Sqrt(μ * *semilatus_rectum) / Radian;
      } else if (specific_angular_momentum) {
        SpecificAngularMomentum const& h = *specific_angular_momentum;
        semilatus_rectum = Pow<2>(h * Radian) / μ;
      }
    }
  };

  complete_conic_parameters(/*first_pass=*/true);

  // TODO(egg): some of these formulae are very ill-conditioned near the
  // parabolic case, and can be easily rewritten.  Investigate.
  if (eccentricity_specified && semimajor_axis_specified) {
    double const& e = *eccentricity;
    Length const& a = *semimajor_axis;
    semiminor_axis = a * Sqrt(1 - Pow<2>(e));
    impact_parameter = -a * Sqrt(Pow<2>(e) - 1);
    semilatus_rectum = a * (1 - Pow<2>(e));
    periapsis_distance = a * (1 - e);
    apoapsis_distance = a * (1 + e);
  }
  if (eccentricity_specified && semiminor_axis_specified) {
    double const& e = *eccentricity;
    Length const& abs_b = e > 1 ? *impact_parameter : *semiminor_axis;
    semilatus_rectum = abs_b * Sqrt(Abs(1 - Pow<2>(e)));
    semimajor_axis = *semilatus_rectum / (Pow<2>(e) - 1);
    periapsis_distance = *semimajor_axis * (1 - e);
    apoapsis_distance = *semimajor_axis * (1 + e);
  }
  if (eccentricity_specified && semilatus_rectum_specified) {
    double const& e = *eccentricity;
    Length const& p = *semilatus_rectum;
    semimajor_axis = p / (1 - Pow<2>(e));
    semiminor_axis = *semimajor_axis * Sqrt(1 - Pow<2>(e));
    impact_parameter = *semimajor_axis * Sqrt(Pow<2>(e) - 1);
    periapsis_distance = p / (1 + e);
    apoapsis_distance = p / (1 - e);
  }
  if (eccentricity_specified && periapsis_distance_specified) {
    double const& e = *eccentricity;
    Length const& r_pe = *periapsis_distance;
    semimajor_axis = r_pe / (1 - e);
    semiminor_axis = *semimajor_axis * Sqrt(1 - Pow<2>(e));
    impact_parameter = *semimajor_axis * Sqrt(Pow<2>(e) - 1);
    semilatus_rectum = r_pe * (1 + e);
    apoapsis_distance = *semimajor_axis * (1 - e);
  }
  if (eccentricity_specified && apoapsis_distance_specified) {
    double const& e = *eccentricity;
    Length const& r_ap = *apoapsis_distance;
    semimajor_axis = r_ap / (1 + e);
    semiminor_axis = *semimajor_axis * Sqrt(1 - Pow<2>(e));
    impact_parameter = *semimajor_axis * Sqrt(Pow<2>(e) - 1);
    semilatus_rectum = r_ap * (1 - e);
    periapsis_distance = *semimajor_axis * (1 - e);
  }
  if (semimajor_axis_specified && semiminor_axis_specified) {
    Length const& a = *semimajor_axis;
    auto const& b² = *semiminor_axis != *semiminor_axis
                           ? -Pow<2>(*impact_parameter)
                           : Pow<2>(*semiminor_axis);
    eccentricity = Sqrt(1 - b² / Pow<2>(a));
    semilatus_rectum = b² / a;
    periapsis_distance = a - Sqrt(Pow<2>(a) - b²);
    apoapsis_distance = a + Sqrt(Pow<2>(a) - b²);
  }
  if (semimajor_axis_specified && semilatus_rectum_specified) {
    Length const& a = *semimajor_axis;
    Length const& p = *semilatus_rectum;
    eccentricity = Sqrt(1 - p / a);
    semiminor_axis = Sqrt(a * p);
    impact_parameter = Sqrt(-a * p);
    double const& e = *eccentricity;
    periapsis_distance = p / (1 + e);
    apoapsis_distance = p / (1 - e);
  }
  if (semimajor_axis_specified && periapsis_distance_specified) {
    Length const& a = *semimajor_axis;
    Length const& r_pe = *periapsis_distance;
    eccentricity = 1 - r_pe / a;
    semiminor_axis = Sqrt((2 * a - r_pe) * r_pe);
    impact_parameter = Sqrt((r_pe - 2 * a) * r_pe);
    semilatus_rectum = r_pe * (2 - r_pe / a);
    apoapsis_distance = 2 * a - r_pe;
  }
  if (semimajor_axis_specified && apoapsis_distance_specified) {
    Length const& a = *semimajor_axis;
    Length const& r_ap = *apoapsis_distance;
    eccentricity = r_ap / a - 1;
    semiminor_axis = Sqrt((2 * a - r_ap) * r_ap);
    impact_parameter = Sqrt((r_ap - 2 * a) * r_ap);
    semilatus_rectum = r_ap * (2 - r_ap / a);
    periapsis_distance = 2 * a - r_ap;
  }
  if (semiminor_axis_specified && semilatus_rectum_specified) {
    auto const& b² = *semiminor_axis != *semiminor_axis
                           ? -Pow<2>(*impact_parameter)
                           : Pow<2>(*semiminor_axis);
    Length const& p = *semilatus_rectum;
    eccentricity = Sqrt(1 - Pow<2>(p) / b²);
    semimajor_axis = b² / p;
    double const& e = *eccentricity;
    periapsis_distance = p / (1 + e);
    apoapsis_distance = p / (1 - e);
  }
  if (semiminor_axis_specified && periapsis_distance_specified) {
    auto const& b² = *semiminor_axis != *semiminor_axis
                           ? -Pow<2>(*impact_parameter)
                           : Pow<2>(*semiminor_axis);
    Length const& r_pe = *periapsis_distance;
    eccentricity = (b² - Pow<2>(r_pe)) / (b² + Pow<2>(r_pe));
    semimajor_axis = (b² + Pow<2>(r_pe)) / (2 * r_pe);
    semilatus_rectum = 2 * b² * r_pe / (b² + Pow<2>(r_pe));
    apoapsis_distance = b² / r_pe;
  }
  if (semiminor_axis_specified && apoapsis_distance_specified) {
    auto const& b² = *semiminor_axis != *semiminor_axis
                           ? -Pow<2>(*impact_parameter)
                           : Pow<2>(*semiminor_axis);
    Length const& r_ap = *apoapsis_distance;
    eccentricity = (Pow<2>(r_ap) - b²) / (b² + Pow<2>(r_ap));
    semimajor_axis = (b² + Pow<2>(r_ap)) / (2 * r_ap);
    semilatus_rectum = 2 * b² * r_ap / (b² + Pow<2>(r_ap));
    periapsis_distance = b² / r_ap;
  }
  if (semilatus_rectum_specified && periapsis_distance_specified) {
    Length const& p = *semilatus_rectum;
    Length const& r_pe = *periapsis_distance;
    eccentricity = p / r_pe - 1;
    semimajor_axis = Pow<2>(r_pe) / (2 * r_pe - p);
    Length const& a = *semimajor_axis;
    semiminor_axis = Sqrt(p * a);
    impact_parameter = Sqrt(-p * a);
    apoapsis_distance = p * r_pe / (2 * r_pe - p);
  }
  if (semilatus_rectum_specified && apoapsis_distance_specified) {
    Length const& p = *semilatus_rectum;
    Length const& r_ap = *apoapsis_distance;
    eccentricity = 1 - p / r_ap;
    semimajor_axis = Pow<2>(r_ap) / (2 * r_ap - p);
    Length const& a = *semimajor_axis;
    semiminor_axis = Sqrt(p * a);
    impact_parameter = Sqrt(-p * a);
    periapsis_distance = p * r_ap / (2 * r_ap - p);
  }
  if (periapsis_distance_specified && periapsis_distance_specified) {
    Length const& r_pe = *periapsis_distance;
    Length const& r_ap = *apoapsis_distance;
    eccentricity = (r_ap - r_pe) / (r_ap + r_pe);
    semimajor_axis = (r_ap + r_pe) / 2;
    semiminor_axis = Sqrt(r_ap * r_pe);
    impact_parameter = Sqrt(-r_ap * r_pe);
    semilatus_rectum = (2 * r_ap * r_pe) / (r_ap + r_pe);
  }

  complete_conic_parameters(/*first_pass=*/false);

  auto& argument_of_periapsis = elements.argument_of_periapsis;
  auto& longitude_of_periapsis = elements.longitude_of_periapsis;
  auto const& Ω = elements.longitude_of_ascending_node;

  CHECK_EQ(
      argument_of_periapsis.has_value() + longitude_of_periapsis.has_value(),
      1);

  if (argument_of_periapsis) {
    auto const& ω = *argument_of_periapsis;
    longitude_of_periapsis = Ω + ω;
  } else if (longitude_of_periapsis) {
    auto const& ϖ = *longitude_of_periapsis;
    argument_of_periapsis = ϖ - Ω;
  }

  CompleteAnomalies(elements);
}

template<typename Frame>
void KeplerOrbit<Frame>::CompleteAnomalies(KeplerianElements<Frame>& elements) {
  auto const& e = *elements.eccentricity;
  auto& true_anomaly = elements.true_anomaly;
  auto& mean_anomaly = elements.mean_anomaly;
  auto& hyperbolic_mean_anomaly = elements.hyperbolic_mean_anomaly;
  CHECK_EQ(true_anomaly.has_value() + mean_anomaly.has_value() +
               hyperbolic_mean_anomaly.has_value(),
           1);
  if (true_anomaly) {
    auto const& ν = *true_anomaly;
    auto const positive_angle = [](Angle const& α) -> Angle {
      return α > 0 * Radian ? α : α + 2 * π * Radian;
    };
    Angle const eccentric_anomaly =
        ArcTan(Sqrt(1 - Pow<2>(e)) * Sin(ν), e + Cos(ν));
    mean_anomaly =
        positive_angle(eccentric_anomaly - e * Sin(eccentric_anomaly) * Radian);
    Angle const hyperbolic_eccentric_anomaly =
        ArcCosh((Cos(ν) - e) / (1 - e * Cos(ν)));
    hyperbolic_mean_anomaly =
        positive_angle(e * Sinh(hyperbolic_eccentric_anomaly) * Radian -
                       hyperbolic_eccentric_anomaly);
  } else if (mean_anomaly) {
    auto const kepler_equation =
        [e, mean_anomaly](Angle const& eccentric_anomaly) -> Angle {
      return *mean_anomaly -
             (eccentric_anomaly - e * Sin(eccentric_anomaly) * Radian);
    };
    Angle const eccentric_anomaly = e == 0 ? *mean_anomaly
                                           : Bisect(kepler_equation,
                                                    *mean_anomaly - e * Radian,
                                                    *mean_anomaly + e * Radian);
    true_anomaly = 2 * ArcTan(Sqrt(1 + e) * Sin(eccentric_anomaly / 2),
                              Sqrt(1 - e) * Cos(eccentric_anomaly / 2));
    hyperbolic_mean_anomaly = NaN<Angle>();
  } else if (hyperbolic_mean_anomaly) {
    auto const hyperbolic_kepler_equation =
        [e, hyperbolic_mean_anomaly](
            Angle const& hyperbolic_eccentric_anomaly) -> Angle {
      return *hyperbolic_mean_anomaly -
             (e * Sinh(hyperbolic_eccentric_anomaly) * Radian -
              hyperbolic_eccentric_anomaly);
    };
    Angle const hyperbolic_eccentric_anomaly =
        Bisect(hyperbolic_kepler_equation,
               *hyperbolic_mean_anomaly - e * Radian,
               *hyperbolic_mean_anomaly + e * Radian);
    true_anomaly =
        2 * ArcTan(Sqrt(1 + e) * Sinh(hyperbolic_eccentric_anomaly / 2),
                   Sqrt(1 - e) * Cosh(hyperbolic_eccentric_anomaly / 2));
    mean_anomaly = NaN<Angle>();
  }
}

}  // namespace internal_kepler_orbit
}  // namespace physics
}  // namespace principia
