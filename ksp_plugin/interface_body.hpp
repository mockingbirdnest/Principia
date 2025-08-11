#pragma once

#include "ksp_plugin/interface.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/strings/str_split.h"
#include "base/array.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "geometry/space_transformations.hpp"
#include "integrators/integrators.hpp"
#include "ksp_plugin/orbit_analyser.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/renderer.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::base::_array;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_space_transformations;
using namespace principia::integrators::_integrators;
using namespace principia::ksp_plugin::_orbit_analyser;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rigid_motion;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

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

template<typename Frame>
struct XYZConverter<AngularVelocity<Frame>> {
  static AngularVelocity<Frame> FromXYZ(XYZ const& xyz) {
    return AngularVelocity<Frame>(interface::FromXYZ(xyz) * (Radian / Second));
  }
  static XYZ ToXYZ(AngularVelocity<Frame> const& velocity) {
    return interface::ToXYZ(velocity.coordinates() / (Radian / Second));
  }
};

template<typename Frame>
struct XYZConverter<Bivector<AngularMomentum, Frame>> {
  static constexpr AngularMomentum mts_unit =
      Pow<2>(Metre) * Tonne * Radian / Second;
  static Bivector<AngularMomentum, Frame> FromXYZ(XYZ const& xyz) {
    return Bivector<AngularMomentum, Frame>(interface::FromXYZ(xyz) * mts_unit);
  }
  static XYZ ToXYZ(Bivector<AngularMomentum, Frame> const& velocity) {
    return interface::ToXYZ(velocity.coordinates() / mts_unit);
  }
};

template<>
struct XYZConverter<R3Element<MomentOfInertia>> {
  static constexpr MomentOfInertia mts_unit = Pow<2>(Metre) * Tonne;
  static R3Element<MomentOfInertia> FromXYZ(XYZ const& xyz) {
    return R3Element<MomentOfInertia>(
        xyz.x * mts_unit, xyz.y * mts_unit, xyz.z * mts_unit);
  }
  static XYZ ToXYZ(R3Element<MomentOfInertia> const& moments_of_inertia) {
    return {moments_of_inertia.x / mts_unit,
            moments_of_inertia.y / mts_unit,
            moments_of_inertia.z / mts_unit};
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

inline bool operator==(Node const& left, Node const& right) {
  return NaNIndependentEq(left.time, right.time) &&
         left.world_position == right.world_position &&
         NaNIndependentEq(left.apparent_inclination_in_degrees,
                          right.apparent_inclination_in_degrees) &&
         NaNIndependentEq(left.out_of_plane_velocity,
                          right.out_of_plane_velocity);
}

inline bool operator==(OrbitAnalysis const& left, OrbitAnalysis const& right) {
  return left.elements == right.elements &&
         left.ground_track_equatorial_crossings ==
             right.ground_track_equatorial_crossings &&
         left.solar_times_of_nodes == right.solar_times_of_nodes &&
         left.mission_duration == right.mission_duration &&
         left.primary_index == right.primary_index &&
         left.progress_of_next_analysis == right.progress_of_next_analysis &&
         left.recurrence == right.recurrence;
}

inline bool operator==(EquatorialCrossings const& left,
                       EquatorialCrossings const& right) {
  return left.longitudes_reduced_to_ascending_pass ==
             right.longitudes_reduced_to_ascending_pass &&
         left.longitudes_reduced_to_descending_pass ==
             right.longitudes_reduced_to_descending_pass;
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

inline bool operator==(QPRW const& left, QPRW const& right) {
  return left.qp == right.qp && left.r == right.r && left.w == right.w;
}

inline bool operator==(SolarTimesOfNodes const& left,
                       SolarTimesOfNodes const& right) {
  return left.mean_solar_times_of_ascending_nodes ==
             right.mean_solar_times_of_ascending_nodes &&
         left.mean_solar_times_of_descending_nodes ==
             right.mean_solar_times_of_descending_nodes;
}

inline bool operator==(Status const& left, Status const& right) {
  return left.error == right.error;
}

inline bool operator==(TQP const& left, TQP const& right) {
  return NaNIndependentEq(left.t, right.t) && left.qp == right.qp;
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

inline physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters
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
    physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters,
    physics::_ephemeris::Ephemeris<
        Barycentric>::GeneralizedAdaptiveStepParameters>
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

inline Renderer::Node FromNode(Plugin const& plugin,
                               Node const& node) {
  return Renderer::Node{
      .time = FromGameTime(plugin, node.time),
      .position = FromXYZ<Position<World>>(node.world_position),
      .apparent_inclination = node.apparent_inclination_in_degrees * Degree,
      .out_of_plane_velocity = node.out_of_plane_velocity * Metre / Second};
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

inline Quaternion FromWXYZ(WXYZ const& wxyz) {
  // It is critical to normalize the quaternion that we receive from Unity: it
  // is normalized in *single* precision, which is fine for KSP where the moving
  // origin of World ensures that coordinates are never very large.  But in the
  // C++ code we do some computations in Barycentric, which typically results in
  // large coordinates for which we need normalization in *double* precision.
  return Normalize(Quaternion{wxyz.w, {wxyz.x, wxyz.y, wxyz.z}});
}

inline R3Element<double> FromXYZ(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

template<>
inline Position<World> FromXYZ<Position<World>>(XYZ const& xyz) {
  return XYZConverter<Position<World>>::FromXYZ(xyz);
}

template<>
inline Position<EccentricPart> FromXYZ<Position<EccentricPart>>(
    XYZ const& xyz) {
  return XYZConverter<Position<EccentricPart>>::FromXYZ(xyz);
}

template<>
Velocity<Frenet<NavigationFrame>>
inline FromXYZ<Velocity<Frenet<NavigationFrame>>>(XYZ const& xyz) {
  return XYZConverter<Velocity<Frenet<NavigationFrame>>>::FromXYZ(xyz);
}

template<>
AngularVelocity<World>
inline FromXYZ<AngularVelocity<World>>(XYZ const& xyz) {
  return XYZConverter<AngularVelocity<World>>::FromXYZ(xyz);
}

template<>
R3Element<MomentOfInertia>
inline FromXYZ<R3Element<MomentOfInertia>>(XYZ const& xyz) {
  return XYZConverter<R3Element<MomentOfInertia>>::FromXYZ(xyz);
}

inline AdaptiveStepParameters ToAdaptiveStepParameters(
    physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters const&
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
    physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    physics::_ephemeris::Ephemeris<
        Barycentric>::GeneralizedAdaptiveStepParameters const&
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
    physics::_kepler_orbit::KeplerianElements<Barycentric> const&
        keplerian_elements) {
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

inline Node ToNode(Plugin const& plugin,
                   Renderer::Node const& node) {
  return Node{
      .time = ToGameTime(plugin, node.time),
      .world_position = ToXYZ(node.position),
      .apparent_inclination_in_degrees = node.apparent_inclination / Degree,
      .out_of_plane_velocity = node.out_of_plane_velocity / (Metre / Second)};
}

inline QP ToQP(DegreesOfFreedom<World> const& dof) {
  return QPConverter<DegreesOfFreedom<World>>::ToQP(dof);
}

inline QP ToQP(RelativeDegreesOfFreedom<AliceSun> const& relative_dof) {
  return QPConverter<RelativeDegreesOfFreedom<AliceSun>>::ToQP(relative_dof);
}

inline Status* ToNewStatus(absl::Status const& status) {
  if (status.ok()) {
    return new Status{static_cast<int>(status.code()),
                      /*message=*/nullptr};
  } else {
    std::string_view const message = status.message();
    LOG(ERROR) << message;
    UniqueArray<char> allocated_message(message.size() + 1);
    std::memcpy(allocated_message.data.get(),
                message.data(),
                message.size() + 1);
    return new Status{static_cast<int>(status.code()),
                      allocated_message.data.release()};
  }
}

inline WXYZ ToWXYZ(Quaternion const& quaternion) {
  return {quaternion.real_part(),
          quaternion.imaginary_part().x,
          quaternion.imaginary_part().y,
          quaternion.imaginary_part().z};
}

inline XY ToXY(RP2Point<Length, Camera> const& rp2_point) {
  return {rp2_point.x() / Metre, rp2_point.y() / Metre};
}

inline XYZ ToXYZ(R3Element<double> const& r3_element) {
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

inline XYZ ToXYZ(Bivector<AngularMomentum, World> const& angular_momentum) {
  return XYZConverter<Bivector<AngularMomentum, World>>::ToXYZ(
      angular_momentum);
}

template<typename T>
Interval ToInterval(geometry::_interval::Interval<T> const& interval) {
  return {interval.min / si::Unit<T>,
          interval.max / si::Unit<T>};
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
    case serialization::BarycentricRotatingReferenceFrame::
        kExtensionFieldNumber:
      return plugin.NewBarycentricRotatingNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredBodyDirectionReferenceFrame::
        kExtensionFieldNumber:
      return plugin.NewBodyCentredBodyDirectionNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredNonRotatingReferenceFrame::
        kExtensionFieldNumber:
      return plugin.NewBodyCentredNonRotatingNavigationFrame(
          parameters.centre_index);
    case serialization::BodySurfaceReferenceFrame::kExtensionFieldNumber:
      return plugin.NewBodySurfaceNavigationFrame(parameters.centre_index);
    default:
      LOG(FATAL) << "Unexpected extension " << parameters.extension;
      std::abort();
  }
}

inline not_null<std::unique_ptr<PlottingFrame>> NewPlottingFrame(
    Plugin const& plugin,
    PlottingFrameParameters const& parameters) {
  CHECK_NOTNULL(parameters.primary_index);
  CHECK_NOTNULL(parameters.secondary_index);
  switch (parameters.extension) {
    case serialization::RotatingPulsatingReferenceFrame::
        kExtensionFieldNumber: {
      std::vector<int> primary_indices;
      for (int const* const* index_ptr = parameters.primary_index;
           *index_ptr != nullptr;
           ++index_ptr) {
        primary_indices.push_back(**index_ptr);
      }
      std::vector<int> secondary_indices;
      for (int const* const* index_ptr = parameters.secondary_index;
           *index_ptr != nullptr;
           ++index_ptr) {
        secondary_indices.push_back(**index_ptr);
      }
      return plugin.NewRotatingPulsatingPlottingFrame(primary_indices,
                                                      secondary_indices);
    }
    default:
      int primary_index;
      if (*parameters.primary_index == nullptr) {
        primary_index = -1;
      } else {
        primary_index = **parameters.primary_index;
      }
      int secondary_index;
      if (*parameters.secondary_index == nullptr) {
        secondary_index = -1;
      } else {
        secondary_index = **parameters.secondary_index;
      }
      return NewNavigationFrame(plugin,
                                {.extension = parameters.extension,
                                 .centre_index = parameters.centre_index,
                                 .primary_index = primary_index,
                                 .secondary_index = secondary_index});
  }
}

inline RigidMotion<EccentricPart, World> MakePartRigidMotion(
    QP const& part_world_degrees_of_freedom,
    WXYZ const& part_rotation,
    XYZ const& part_angular_velocity) {
  DegreesOfFreedom<World> const part_degrees_of_freedom =
      FromQP<DegreesOfFreedom<World>>(part_world_degrees_of_freedom);
  Rotation<EccentricPart, World> const part_to_world(FromWXYZ(part_rotation));
  RigidTransformation<EccentricPart, World> const part_rigid_transformation(
      EccentricPart::origin,
      part_degrees_of_freedom.position(),
      part_to_world.Forget<OrthogonalMap>());
  RigidMotion<EccentricPart, World> part_rigid_motion(
      part_rigid_transformation,
      FromXYZ<AngularVelocity<World>>(part_angular_velocity),
      part_degrees_of_freedom.velocity());
  return part_rigid_motion;
}

// Same as `MakePartRigidMotion`, but uses the separate type `ApparentWorld` to
// avoid mixing uncorrected and corrected data.
inline RigidMotion<EccentricPart, ApparentWorld> MakePartApparentRigidMotion(
    QP const& part_world_degrees_of_freedom,
    WXYZ const& part_rotation,
    XYZ const& part_angular_velocity) {
  return RigidMotion<World, ApparentWorld>::Identity() *
         MakePartRigidMotion(part_world_degrees_of_freedom,
                             part_rotation,
                             part_angular_velocity);
}

// Ownership is returned to the caller.  Note that the result may own additional
// objects via C pointers; it must be not be deleted from C++, and must instead
// be passed to the generated C# marshaller, which will properly delete it.
inline not_null<OrbitAnalysis*> NewOrbitAnalysis(
    OrbitAnalyser::Analysis* const vessel_analysis,
    Plugin const& plugin,
    int const* const revolutions_per_cycle,
    int const* const days_per_cycle,
    int const ground_track_revolution) {
  auto* const analysis = new OrbitAnalysis{};
  CHECK_EQ(revolutions_per_cycle == nullptr, days_per_cycle == nullptr);
  bool const has_nominal_recurrence = revolutions_per_cycle != nullptr;
  if (has_nominal_recurrence) {
    CHECK_GT(*revolutions_per_cycle, 0);
    CHECK_NE(*days_per_cycle, 0);
  }
  if (vessel_analysis == nullptr) {
    return analysis;
  }
  analysis->primary_index =
      vessel_analysis->primary() == nullptr
          ? nullptr
          : new int(plugin.CelestialIndexOfBody(*vessel_analysis->primary()));

  analysis->mission_duration = vessel_analysis->mission_duration() / Second;
  if (vessel_analysis->elements().has_value()) {
    auto const& elements = *vessel_analysis->elements();
    analysis->elements = new OrbitalElements{
        .sidereal_period = elements.sidereal_period() / Second,
        .nodal_period = elements.nodal_period() / Second,
        .anomalistic_period = elements.anomalistic_period() / Second,
        .nodal_precession = elements.nodal_precession() / (Radian / Second),
        .mean_semimajor_axis =
            ToInterval(elements.mean_semimajor_axis_interval()),
        .mean_eccentricity = ToInterval(elements.mean_eccentricity_interval()),
        .mean_inclination = ToInterval(elements.mean_inclination_interval()),
        .mean_longitude_of_ascending_nodes =
            ToInterval(elements.mean_longitude_of_ascending_node_interval()),
        .mean_argument_of_periapsis =
            ToInterval(elements.mean_argument_of_periapsis_interval()),
        .mean_periapsis_distance =
            ToInterval(elements.mean_periapsis_distance_interval()),
        .mean_apoapsis_distance =
            ToInterval(elements.mean_apoapsis_distance_interval()),
        .radial_distance =
            ToInterval(*vessel_analysis->radial_distance_interval()),
    };
  }
  if (has_nominal_recurrence && vessel_analysis->primary() != nullptr) {
    int const Cᴛₒ =
        Sign(vessel_analysis->primary()->angular_frequency()) *
        std::abs(*days_per_cycle);
    int const νₒ =
        std::nearbyint(static_cast<double>(*revolutions_per_cycle) / Cᴛₒ);
    int const Dᴛₒ = *revolutions_per_cycle - νₒ * Cᴛₒ;
    int const gcd = std::gcd(Dᴛₒ, Cᴛₒ);
    vessel_analysis->SetRecurrence({νₒ, Dᴛₒ / gcd, Cᴛₒ / gcd});
  } else {
    vessel_analysis->ResetRecurrence();
  }
  if (vessel_analysis->recurrence().has_value()) {
    auto const& recurrence = *vessel_analysis->recurrence();
    analysis->recurrence = new OrbitRecurrence{
        .nuo = recurrence.νₒ(),
        .dto = recurrence.Dᴛₒ(),
        .cto = recurrence.Cᴛₒ(),
        .number_of_revolutions = recurrence.number_of_revolutions(),
        .equatorial_shift = recurrence.equatorial_shift() / Radian,
        .base_interval = recurrence.base_interval() / Radian,
        .grid_interval = recurrence.grid_interval() / Radian,
        .subcycle = recurrence.subcycle(),
    };
  }
  if (auto const& ground_track = vessel_analysis->ground_track();
      ground_track.has_value()) {
    if (ground_track->mean_solar_times_of_ascending_nodes().has_value() &&
        ground_track->mean_solar_times_of_descending_nodes().has_value()) {
      analysis->solar_times_of_nodes = new SolarTimesOfNodes{
          .mean_solar_times_of_ascending_nodes =
              ToInterval(*ground_track->mean_solar_times_of_ascending_nodes()),
          .mean_solar_times_of_descending_nodes = ToInterval(
              *ground_track->mean_solar_times_of_descending_nodes())};
    }
    if (vessel_analysis->equatorial_crossings().has_value()) {
      auto const& equatorial_crossings =
          *vessel_analysis->equatorial_crossings();
      analysis->ground_track_equatorial_crossings = new EquatorialCrossings{
          .longitudes_reduced_to_ascending_pass =
              ToInterval(equatorial_crossings.longitudes_reduced_to_pass(
                  2 * ground_track_revolution - 1)),
          .longitudes_reduced_to_descending_pass =
              ToInterval(equatorial_crossings.longitudes_reduced_to_pass(
                  2 * ground_track_revolution)),
      };
    }
  }
  return analysis;
}

inline FlightPlan& GetFlightPlan(Plugin const& plugin,
                                 char const* const vessel_guid) {
  Vessel& vessel = *plugin.GetVessel(vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  // Force deserialization of the flight plan, now that we actually need it.
  vessel.ReadFlightPlanFromMessage();
  return vessel.flight_plan();
}

}  // namespace interface
}  // namespace principia
