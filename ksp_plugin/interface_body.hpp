
#pragma once

#include "ksp_plugin/interface.hpp"

#include <cmath>
#include <limits>
#include <utility>

namespace principia {
namespace interface {

using integrators::AdaptiveStepSizeIntegrator;
using physics::Ephemeris;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

// No partial specialization of functions, so we wrap everything into structs.
// C++, I hate you.

template<typename T>
struct QPConverter {};

template<typename T>
struct XYZConverter {};

template<typename Frame>
struct QPConverter<DegreesOfFreedom<Frame>> {
  static DegreesOfFreedom<Frame> FromQP(QP const& qp) {
    return DegreesOfFreedom<Frame>(
        XYZConverter<Position<Frame>>::FromXYZ(qp.q),
        XYZConverter<Velocity<Frame>>::FromXYZ(qp.p));
  }
  static QP ToQP(DegreesOfFreedom<Frame> const& dof) {
    return {XYZConverter<Position<Frame>>::ToXYZ(dof.position()),
            XYZConverter<Velocity<Frame>>::ToXYZ(dof.velocity())};
  }
};

template<typename Frame>
struct QPConverter<RelativeDegreesOfFreedom<Frame>> {
  static RelativeDegreesOfFreedom<Frame> FromQP(QP const& qp) {
    return RelativeDegreesOfFreedom<Frame>(
        XYZConverter<Displacement<Frame>>::FromXYZ(qp.q),
        XYZConverter<Velocity<Frame>>::FromXYZ(qp.p));
  }
  static QP ToQP(RelativeDegreesOfFreedom<Frame> const& relative_dof) {
    return {
        XYZConverter<Displacement<Frame>>::ToXYZ(relative_dof.displacement()),
        XYZConverter<Velocity<Frame>>::ToXYZ(relative_dof.velocity())};
  }
};

template<typename Frame>
struct XYZConverter<Displacement<Frame>> {
  static Displacement<Frame> FromXYZ(XYZ const& xyz) {
    return Displacement<Frame>(interface::FromXYZ(xyz) * Metre);
  }
  static XYZ ToXYZ(Displacement<Frame> const& displacement) {
    return interface::ToXYZ(displacement.coordinates() / Metre);
  }
};

template<typename Frame>
struct XYZConverter<Position<Frame>> {
  static Position<Frame> FromXYZ(XYZ const& xyz) {
    return Position<Frame>(Frame::origin +
                           XYZConverter<Displacement<Frame>>::FromXYZ(xyz));
  }
  static XYZ ToXYZ(Position<Frame> const& position) {
    return XYZConverter<Displacement<Frame>>::ToXYZ(position - Frame::origin);
  }
};

template<typename Frame>
struct XYZConverter<Velocity<Frame>> {
  static Velocity<Frame> FromXYZ(XYZ const& xyz) {
    return Velocity<Frame>(interface::FromXYZ(xyz) * (Metre / Second));
  }
  static XYZ ToXYZ(Velocity<Frame> const& velocity) {
    return interface::ToXYZ(velocity.coordinates() / (Metre / Second));
  }
};

inline bool NaNIndependentEq(double const left, double const right) {
  return (left == right) || (std::isnan(left) && std::isnan(right));
}

template<typename T>
std::unique_ptr<T> TakeOwnership(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

template<typename T>
std::unique_ptr<T[]> TakeOwnershipArray(T** const pointer) {
  CHECK_NOTNULL(pointer);
  std::unique_ptr<T[]> owned_pointer(*pointer);
  *pointer = nullptr;
  return owned_pointer;
}

inline bool operator==(AdaptiveStepParameters const& left,
                       AdaptiveStepParameters const& right) {
  return left.integrator_kind == right.integrator_kind &&
         left.max_steps == right.max_steps &&
         NaNIndependentEq(left.length_integration_tolerance,
                          right.length_integration_tolerance) &&
         NaNIndependentEq(left.speed_integration_tolerance,
                          right.speed_integration_tolerance);
}

inline bool operator==(Burn const& left, Burn const& right) {
  return NaNIndependentEq(left.thrust_in_kilonewtons,
                          right.thrust_in_kilonewtons) &&
         NaNIndependentEq(left.specific_impulse_in_seconds_g0,
                          right.specific_impulse_in_seconds_g0) &&
         left.frame == right.frame &&
         NaNIndependentEq(left.initial_time, right.initial_time) &&
         left.delta_v == right.delta_v;
}

inline bool operator==(FlightPlanAdaptiveStepParameters const& left,
                       FlightPlanAdaptiveStepParameters const& right) {
  return left.integrator_kind == right.integrator_kind &&
         left.generalized_integrator_kind ==
             right.generalized_integrator_kind &&
         left.max_steps == right.max_steps &&
         NaNIndependentEq(left.length_integration_tolerance,
                          right.length_integration_tolerance) &&
         NaNIndependentEq(left.speed_integration_tolerance,
                          right.speed_integration_tolerance);
}

inline bool operator==(Interval const& left, Interval const& right) {
  return NaNIndependentEq(left.min, right.min) &&
         NaNIndependentEq(left.max, right.max);
}

inline bool operator==(NavigationFrameParameters const& left,
                       NavigationFrameParameters const& right) {
  return left.extension == right.extension &&
         left.centre_index == right.centre_index &&
         left.primary_index == right.primary_index &&
         left.secondary_index == right.secondary_index;
}

inline bool operator==(NavigationManoeuvre const& left,
                       NavigationManoeuvre const& right) {
  return left.burn == right.burn &&
         NaNIndependentEq(left.initial_mass_in_tonnes,
                          right.initial_mass_in_tonnes) &&
         NaNIndependentEq(left.final_mass_in_tonnes,
                          right.final_mass_in_tonnes) &&
         NaNIndependentEq(left.mass_flow, right.mass_flow) &&
         NaNIndependentEq(left.duration, right.duration) &&
         NaNIndependentEq(left.final_time, right.final_time) &&
         NaNIndependentEq(left.time_of_half_delta_v,
                          right.time_of_half_delta_v) &&
         NaNIndependentEq(left.time_to_half_delta_v,
                          right.time_to_half_delta_v);
}

inline bool operator==(NavigationManoeuvreFrenetTrihedron const& left,
                       NavigationManoeuvreFrenetTrihedron const& right) {
  return left.binormal == right.binormal &&
         left.normal == right.normal &&
         left.tangent == right.tangent;
}

inline bool operator==(OrbitAnalysis const& left, OrbitAnalysis const& right) {
  return left.elements == right.elements &&
         left.ground_track == right.ground_track &&
         left.mission_duration == right.mission_duration &&
         left.primary_index == right.primary_index &&
         left.progress_percentage && right.progress_percentage &&
         left.recurrence == right.recurrence;
}

inline bool operator==(EquatorialCrossings const& left,
                       EquatorialCrossings const& right) {
  return left.longitudes_reduced_to_ascending_pass ==
             right.longitudes_reduced_to_ascending_pass &&
         left.longitudes_reduced_to_descending_pass ==
             right.longitudes_reduced_to_descending_pass;
}

inline bool operator==(OrbitGroundTrack const& left,
                       OrbitGroundTrack const& right) {
  return left.equatorial_crossings ==
             right.equatorial_crossings;
}

inline bool operator==(OrbitRecurrence const& left,
                       OrbitRecurrence const& right) {
  return NaNIndependentEq(left.base_interval, right.base_interval) &&
         left.cto == right.cto && left.dto == right.dto &&
         NaNIndependentEq(left.equatorial_shift, right.equatorial_shift) &&
         NaNIndependentEq(left.grid_interval, right.grid_interval) &&
         left.number_of_revolutions == right.number_of_revolutions &&
         left.nuo == right.nuo && left.subcycle == right.subcycle;
}

inline bool operator==(OrbitalElements const& left,
                       OrbitalElements const& right) {
  return NaNIndependentEq(left.anomalistic_period, right.anomalistic_period) &&
         left.mean_argument_of_periapsis == right.mean_argument_of_periapsis &&
         left.mean_eccentricity == right.mean_eccentricity &&
         left.mean_inclination == right.mean_inclination &&
         left.mean_longitude_of_ascending_nodes ==
             right.mean_longitude_of_ascending_nodes &&
         left.mean_semimajor_axis == right.mean_semimajor_axis &&
         NaNIndependentEq(left.nodal_period, right.nodal_period) &&
         NaNIndependentEq(left.nodal_precession, right.nodal_precession) &&
         NaNIndependentEq(left.sidereal_period, right.sidereal_period);
}

inline bool operator==(QP const& left, QP const& right) {
  return left.q == right.q && left.p == right.p;
}

inline bool operator==(Status const& left, Status const& right) {
  return left.error == right.error;
}

inline bool operator==(WXYZ const& left, WXYZ const& right) {
  return NaNIndependentEq(left.w, right.w) &&
         NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y) &&
         NaNIndependentEq(left.z, right.z);
}

inline bool operator==(XY const& left, XY const& right) {
  return NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y);
}

inline bool operator==(XYZ const& left, XYZ const& right) {
  return NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y) &&
         NaNIndependentEq(left.z, right.z);
}

inline physics::Ephemeris<Barycentric>::AdaptiveStepParameters
FromAdaptiveStepParameters(
    AdaptiveStepParameters const& adaptive_step_parameters) {
  serialization::AdaptiveStepSizeIntegrator message;
  CHECK(serialization::AdaptiveStepSizeIntegrator::Kind_IsValid(
      adaptive_step_parameters.integrator_kind));
  message.set_kind(static_cast<serialization::AdaptiveStepSizeIntegrator::Kind>(
      adaptive_step_parameters.integrator_kind));
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      AdaptiveStepSizeIntegrator<Ephemeris<
          Barycentric>::NewtonianMotionEquation>::ReadFromMessage(message),
      adaptive_step_parameters.max_steps,
      adaptive_step_parameters.length_integration_tolerance * Metre,
      adaptive_step_parameters.speed_integration_tolerance * (Metre / Second));
}

inline std::pair<
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters,
    physics::Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters>
FromFlightPlanAdaptiveStepParameters(FlightPlanAdaptiveStepParameters const&
                                         flight_plan_adaptive_step_parameters) {
  serialization::AdaptiveStepSizeIntegrator message1;
  CHECK(serialization::AdaptiveStepSizeIntegrator::Kind_IsValid(
      flight_plan_adaptive_step_parameters.integrator_kind));
  message1.set_kind(
      static_cast<serialization::AdaptiveStepSizeIntegrator::Kind>(
          flight_plan_adaptive_step_parameters.integrator_kind));

  serialization::AdaptiveStepSizeIntegrator message2;
  CHECK(serialization::AdaptiveStepSizeIntegrator::Kind_IsValid(
      flight_plan_adaptive_step_parameters.generalized_integrator_kind));
  message2.set_kind(
      static_cast<serialization::AdaptiveStepSizeIntegrator::Kind>(
          flight_plan_adaptive_step_parameters.generalized_integrator_kind));

  auto const max_steps = flight_plan_adaptive_step_parameters.max_steps;
  auto const length_integration_tolerance =
      flight_plan_adaptive_step_parameters.length_integration_tolerance * Metre;
  auto const speed_integration_tolerance =
      flight_plan_adaptive_step_parameters.speed_integration_tolerance *
      (Metre / Second);

  auto const adaptive_step_parameters =
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          AdaptiveStepSizeIntegrator<Ephemeris<Barycentric>::
              NewtonianMotionEquation>::ReadFromMessage(message1),
          max_steps,
          length_integration_tolerance,
          speed_integration_tolerance);
  auto const generalized_adaptive_step_parameters =
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
          AdaptiveStepSizeIntegrator<Ephemeris<Barycentric>::
              GeneralizedNewtonianMotionEquation>::ReadFromMessage(message2),
          max_steps,
          length_integration_tolerance,
          speed_integration_tolerance);
  return {adaptive_step_parameters, generalized_adaptive_step_parameters};
}

template<>
inline DegreesOfFreedom<World> FromQP(QP const& qp) {
  return QPConverter<DegreesOfFreedom<World>>::FromQP(qp);
}

template<>
inline RelativeDegreesOfFreedom<AliceSun> FromQP(QP const& qp) {
  return QPConverter<RelativeDegreesOfFreedom<AliceSun>>::FromQP(qp);
}

template<>
inline RelativeDegreesOfFreedom<World> FromQP(QP const& qp) {
  return QPConverter<RelativeDegreesOfFreedom<World>>::FromQP(qp);
}

inline R3Element<double> FromXYZ(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

template<>
inline Position<World> FromXYZ<Position<World>>(XYZ const& xyz) {
  return XYZConverter<Position<World>>::FromXYZ(xyz);
}

template<>
Velocity<Frenet<NavigationFrame>>
inline FromXYZ<Velocity<Frenet<NavigationFrame>>>(XYZ const& xyz) {
  return XYZConverter<Velocity<Frenet<NavigationFrame>>>::FromXYZ(xyz);
}

inline AdaptiveStepParameters ToAdaptiveStepParameters(
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  serialization::AdaptiveStepSizeIntegrator message;
  adaptive_step_parameters.integrator().WriteToMessage(&message);
  return {message.kind(),
          adaptive_step_parameters.max_steps(),
          adaptive_step_parameters.length_integration_tolerance() / Metre,
          adaptive_step_parameters.speed_integration_tolerance() /
              (Metre / Second)};
}

inline FlightPlanAdaptiveStepParameters ToFlightPlanAdaptiveStepParameters(
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    physics::Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        generalized_adaptive_step_parameters) {
  serialization::AdaptiveStepSizeIntegrator message1;
  adaptive_step_parameters.integrator().WriteToMessage(&message1);
  serialization::AdaptiveStepSizeIntegrator message2;
  generalized_adaptive_step_parameters.integrator().WriteToMessage(&message2);
  // TODO(phl): Should we CHECK that the fields are consistent?
  return {message1.kind(),
          message2.kind(),
          adaptive_step_parameters.max_steps(),
          adaptive_step_parameters.length_integration_tolerance() / Metre,
          adaptive_step_parameters.speed_integration_tolerance() /
              (Metre / Second)};
}

inline KeplerianElements ToKeplerianElements(
    physics::KeplerianElements<Barycentric> const& keplerian_elements) {
  return {*keplerian_elements.eccentricity,
          keplerian_elements.semimajor_axis
              ? *keplerian_elements.semimajor_axis / Metre
              : std::numeric_limits<double>::quiet_NaN(),
          keplerian_elements.mean_motion
              ? *keplerian_elements.mean_motion / (Radian / Second)
              : std::numeric_limits<double>::quiet_NaN(),
          keplerian_elements.inclination / Degree,
          keplerian_elements.longitude_of_ascending_node / Degree,
          *keplerian_elements.argument_of_periapsis / Degree,
          *keplerian_elements.mean_anomaly / Radian};
}

inline QP ToQP(DegreesOfFreedom<World> const& dof) {
  return QPConverter<DegreesOfFreedom<World>>::ToQP(dof);
}

inline QP ToQP(RelativeDegreesOfFreedom<AliceSun> const& relative_dof) {
  return QPConverter<RelativeDegreesOfFreedom<AliceSun>>::ToQP(relative_dof);
}

inline Status ToStatus(base::Status const& status) {
  if (!status.ok()) {
    LOG(ERROR) << status.message();
  }
  return {static_cast<int>(status.error())};
}

inline WXYZ ToWXYZ(geometry::Quaternion const& quaternion) {
  return {quaternion.real_part(),
          quaternion.imaginary_part().x,
          quaternion.imaginary_part().y,
          quaternion.imaginary_part().z};
}

inline XY ToXY(geometry::RP2Point<Length, Camera> const& rp2_point) {
  return {rp2_point.x() / Metre, rp2_point.y() / Metre};
}

inline XYZ ToXYZ(geometry::R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
}

inline XYZ ToXYZ(Position<World> const& position) {
  return XYZConverter<Position<World>>::ToXYZ(position);
}

inline XYZ ToXYZ(Velocity<Frenet<NavigationFrame>> const& velocity) {
  return XYZConverter<Velocity<Frenet<NavigationFrame>>>::ToXYZ(velocity);
}

inline XYZ ToXYZ(Vector<double, World> const& direction) {
  return ToXYZ(direction.coordinates());
}

inline XYZ ToXYZ(Velocity<World> const& velocity) {
  return XYZConverter<Velocity<World>>::ToXYZ(velocity);
}

template<typename T>
Interval ToInterval(geometry::Interval<T> const& interval) {
  return {interval.min / quantities::SIUnit<T>(),
          interval.max / quantities::SIUnit<T>()};
}

inline Instant FromGameTime(Plugin const& plugin,
                            double const t) {
  return plugin.GameEpoch() + t * Second;
}

inline double ToGameTime(Plugin const& plugin,
                         Instant const& t) {
  return (t - plugin.GameEpoch()) / Second;
}

inline not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const& plugin,
    NavigationFrameParameters const& parameters) {
  switch (parameters.extension) {
    case serialization::BarycentricRotatingDynamicFrame::
        kExtensionFieldNumber:
      return plugin.NewBarycentricRotatingNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredBodyDirectionDynamicFrame::
        kExtensionFieldNumber:
      return plugin.NewBodyCentredBodyDirectionNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredNonRotatingDynamicFrame::
        kExtensionFieldNumber:
      return plugin.NewBodyCentredNonRotatingNavigationFrame(
          parameters.centre_index);
    case serialization::BodySurfaceDynamicFrame::kExtensionFieldNumber:
      return plugin.NewBodySurfaceNavigationFrame(parameters.centre_index);
    default:
      LOG(FATAL) << "Unexpected extension " << parameters.extension;
      base::noreturn();
  }
}

}  // namespace interface
}  // namespace principia
