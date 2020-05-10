
#include "ksp_plugin/interface.hpp"

#include <limits>
#include <string>

#include "absl/strings/str_cat.h"
#include "base/array.hpp"
#include "base/status.hpp"
#include "base/status_or.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "physics/apsides.hpp"

namespace principia {
namespace interface {

using base::Error;
using base::StatusOr;
using base::UniqueArray;
using geometry::AngularVelocity;
using geometry::Frame;
using geometry::RadiusLatitudeLongitude;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::Vessel;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::ComputeApsides;
using physics::DiscreteTrajectory;
using physics::OblateBody;
using physics::RigidMotion;
using physics::RigidTransformation;

namespace {

// A wrapper for |MakeStatus(base::Status(error, message))|, since we often
// construct the status in the interface itself.
Status MakeStatus(Error const error, std::string const& message) {
  return ToStatus(base::Status(error, message));
}

Status OK() {
  return ToStatus(base::Status::OK);
}

}  // namespace

Status __cdecl principia__ExternalCelestialGetPosition(
    Plugin const* const plugin,
    int const body_index,
    double const time,
    XYZ* const position) {
  journal::Method<journal::ExternalCelestialGetPosition> m{
      {plugin,
       body_index,
       time},
      {position}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  Instant const t = FromGameTime(*plugin, time);
  if (!plugin->HasCelestial(body_index)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        absl::StrCat("No celestial with index ", body_index)));
  }
  auto const& celestial = plugin->GetCelestial(body_index);
  auto const& trajectory = celestial.trajectory();
  if (t < trajectory.t_min() || t > trajectory.t_max()) {
    return m.Return(
        MakeStatus(Error::OUT_OF_RANGE,
                   (std::stringstream{}
                    << "|time| " << t << " does not lie within the domain ["
                    << trajectory.t_min() << ", " << trajectory.t_max()
                    << "] of the trajectory of " << celestial.body()->name())
                       .str()));
  }
  auto const from_solar_system_barycentre =
      plugin->renderer().BarycentricToWorldSun(plugin->PlanetariumRotation())(
          trajectory.EvaluatePosition(t) - Barycentric::origin);
  *position = ToXYZ(from_solar_system_barycentre.coordinates() / Metre);
  return m.Return(OK());
}

Status __cdecl principia__ExternalCelestialGetSurfacePosition(
    Plugin const* const plugin,
    int const body_index,
    double const planetocentric_latitude_in_degrees,
    double const planetocentric_longitude_in_degrees,
    double const radius,
    double const time,
    XYZ* const position) {
  journal::Method<journal::ExternalCelestialGetSurfacePosition> m{
      {plugin,
       body_index,
       planetocentric_latitude_in_degrees,
       planetocentric_longitude_in_degrees,
       radius,
       time},
      {position}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  Instant const t = FromGameTime(*plugin, time);
  if (!plugin->HasCelestial(body_index)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        absl::StrCat("No celestial with index ", body_index)));
  }
  auto const& celestial = plugin->GetCelestial(body_index);
  auto const& trajectory = celestial.trajectory();
  if (t < trajectory.t_min() || t > trajectory.t_max()) {
    return m.Return(
        MakeStatus(Error::OUT_OF_RANGE,
                   (std::stringstream{}
                    << "|time| " << t << " does not lie within the domain ["
                    << trajectory.t_min() << ", " << trajectory.t_max()
                    << "] of the trajectory of " << celestial.body()->name())
                       .str()));
  }
  using Surface = Frame<enum class SurfaceTag>;
  auto const to_world_axes =
      plugin->renderer().BarycentricToWorldSun(plugin->PlanetariumRotation()) *
      celestial.body()->FromSurfaceFrame<Surface>(t).Forget<OrthogonalMap>();
  auto const planetocentric_displacement = Displacement<Surface>(
      RadiusLatitudeLongitude(radius * Metre,
                              planetocentric_latitude_in_degrees * Degree,
                              planetocentric_longitude_in_degrees * Degree)
          .ToCartesian());
  *position =
      ToXYZ(to_world_axes(planetocentric_displacement).coordinates() / Metre);
  return m.Return(OK());
}

Status __cdecl principia__ExternalFlowFreefall(
    Plugin const* const plugin,
    int const central_body_index,
    QP const world_body_centred_initial_degrees_of_freedom,
    double const t_initial,
    double const t_final,
    QP* const world_body_centred_final_degrees_of_freedom) {
  journal::Method<journal::ExternalFlowFreefall> m{
      {plugin,
       central_body_index,
       world_body_centred_initial_degrees_of_freedom,
       t_initial,
       t_final},
      {world_body_centred_final_degrees_of_freedom}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  return m.Return(MakeStatus(Error::UNIMPLEMENTED,
                             "|ExternalFlowFreefall| is not yet implemented"));
}

Status __cdecl principia__ExternalGeopotentialGetCoefficient(
    Plugin const* const plugin,
    int const body_index,
    int const degree,
    int const order,
    XY* const coefficient) {
  journal::Method<journal::ExternalGeopotentialGetCoefficient> m{
      {plugin,
       body_index,
       degree,
       order},
      {coefficient}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  if (!plugin->HasCelestial(body_index)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        absl::StrCat("No celestial with index ", body_index)));
  }
  if (order < 0 || order > degree) {
    return m.Return(MakeStatus(
        Error::INVALID_ARGUMENT,
        absl::StrCat(u8"Expected 0 ≤ order ≤ degree; got degree = ",
                     degree, ", order = ", order)));
  }
  if (degree == 0) {
    *coefficient = {1, 0};
    return m.Return(OK());
  }
  auto const& body = *plugin->GetCelestial(body_index).body();
  if (!body.is_oblate()) {
    *coefficient = {0, 0};
    return m.Return(OK());
  }
  auto const& oblate_body = dynamic_cast<OblateBody<Barycentric> const&>(body);
  if (degree > oblate_body.geopotential_degree()) {
    *coefficient = {0, 0};
    return m.Return(OK());
  }
  *coefficient = {oblate_body.cos()[degree][order],
                  oblate_body.sin()[degree][order]};
  return m.Return(OK());
}

Status __cdecl principia__ExternalGeopotentialGetReferenceRadius(
    Plugin const* const plugin,
    int const body_index,
    double* const reference_radius) {
  journal::Method<journal::ExternalGeopotentialGetReferenceRadius> m{
      {plugin,
       body_index},
      {reference_radius}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  if (!plugin->HasCelestial(body_index)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        absl::StrCat("No celestial with index ", body_index)));
  }
  auto const& body = *plugin->GetCelestial(body_index).body();
  if (!body.is_oblate()) {
    *reference_radius = body.mean_radius() / Metre;
    return m.Return(OK());
  }
  auto const& oblate_body = dynamic_cast<OblateBody<Barycentric> const&>(body);
  *reference_radius = oblate_body.reference_radius() / Metre;
  return m.Return(OK());
}

Status __cdecl principia__ExternalGetNearestPlannedCoastDegreesOfFreedom(
    Plugin const* const plugin,
    int const central_body_index,
    char const* const vessel_guid,
    int const manoeuvre_index,
    XYZ const world_body_centred_reference_position,
    QP* const world_body_centred_nearest_degrees_of_freedom) {
  journal::Method<journal::ExternalGetNearestPlannedCoastDegreesOfFreedom> m{
      {plugin,
       central_body_index,
       vessel_guid,
       manoeuvre_index,
       world_body_centred_reference_position},
      {world_body_centred_nearest_degrees_of_freedom}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  if (manoeuvre_index < 0) {
    return m.Return(MakeStatus(Error::INVALID_ARGUMENT,
                              "Invalid negative |manoeuvre_index|" +
                                  std::to_string(manoeuvre_index)));
  }
  if (!plugin->HasCelestial(central_body_index)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        "No celestial with index " + std::to_string(central_body_index)));
  }
  if (!plugin->HasVessel(vessel_guid)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        "No vessel with GUID " + std::string(vessel_guid)));
  }
  Vessel const& vessel = *plugin->GetVessel(vessel_guid);
  if (!vessel.has_flight_plan()) {
    return m.Return(MakeStatus(
        Error::FAILED_PRECONDITION,
        "Vessel " + vessel.ShortDebugString() + " has no flight plan"));
  }
  FlightPlan const& flight_plan = vessel.flight_plan();
  if (manoeuvre_index >= flight_plan.number_of_manœuvres()) {
    return m.Return(MakeStatus(
        Error::OUT_OF_RANGE,
        "|manoeuvre_index| " + std::to_string(manoeuvre_index) +
            " out of range, vessel " + vessel.ShortDebugString() + " has " +
            std::to_string(flight_plan.number_of_manœuvres()) +
            u8" planned manœuvres"));
  }
  // The index of the coast segment following the desired manœuvre.
  int const segment_index = manoeuvre_index * 2 + 2;
  if (segment_index >= flight_plan.number_of_segments()) {
    return m.Return(MakeStatus(Error::FAILED_PRECONDITION,
                               u8"A singularity occurs within manœuvre " +
                                   std::to_string(manoeuvre_index) + " of " +
                                   vessel.ShortDebugString()));
  }
  DiscreteTrajectory<Barycentric>::Iterator coast_begin;
  DiscreteTrajectory<Barycentric>::Iterator coast_end;
  flight_plan.GetSegment(segment_index, coast_begin, coast_end);
  auto const body_centred_inertial =
      plugin->NewBodyCentredNonRotatingNavigationFrame(central_body_index);
  DiscreteTrajectory<Navigation> coast;
  for (auto it = coast_begin; it != coast_end; ++it) {
    coast.Append(it->time,
                 body_centred_inertial->ToThisFrameAtTime(it->time)(
                     it->degrees_of_freedom));
  }

  Instant const current_time = plugin->CurrentTime();
  // The given |World| position and requested |World| degrees of freedom are
  // body-centred inertial, so |body_centred_inertial| up to an orthogonal map
  // to world coordinates.  Do the conversion directly.
  // NOTE(egg): it is correct to use the orthogonal map at |current_time|,
  // because |body_centred_inertial| does not rotate with respect to
  // |Barycentric|, so the orthogonal map does not depend on time.
  RigidMotion<Navigation, World> to_world_body_centred_inertial(
      RigidTransformation<Navigation, World>(
          Navigation::origin,
          World::origin,
          plugin->renderer().BarycentricToWorld(plugin->PlanetariumRotation()) *
              body_centred_inertial->FromThisFrameAtTime(
                  current_time).orthogonal_map()),
      Navigation::nonrotating,
      Navigation::unmoving);
  auto const from_world_body_centred_inertial =
      to_world_body_centred_inertial.Inverse();
  Position<Navigation> reference_position =
      from_world_body_centred_inertial.rigid_transformation()(
          FromXYZ<Position<World>>(world_body_centred_reference_position));
  DiscreteTrajectory<Navigation> immobile_reference;
  immobile_reference.Append(coast.front().time,
                            {reference_position, Navigation::unmoving});
  if (coast.Size() > 1) {
    immobile_reference.Append(coast.back().time,
                              {reference_position, Navigation::unmoving});
  }
  DiscreteTrajectory<Navigation> apoapsides;
  DiscreteTrajectory<Navigation> periapsides;
  ComputeApsides(/*reference=*/immobile_reference,
                 coast.begin(),
                 coast.end(),
                 /*max_points=*/std::numeric_limits<int>::max(),
                 apoapsides,
                 periapsides);
  if (periapsides.Empty()) {
    bool const begin_is_nearest =
        (coast.front().degrees_of_freedom.position() -
         reference_position).Norm²() <
        (coast.back().degrees_of_freedom.position() -
         reference_position).Norm²();
    *world_body_centred_nearest_degrees_of_freedom =
        ToQP(to_world_body_centred_inertial(
            begin_is_nearest ? coast.front().degrees_of_freedom
                             : coast.back().degrees_of_freedom));
  } else {
    *world_body_centred_nearest_degrees_of_freedom =
        ToQP(to_world_body_centred_inertial(
            periapsides.front().degrees_of_freedom));
  }
  return m.Return(OK());
}

Status __cdecl principia__ExternalVesselGetPosition(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const time,
    XYZ* const position) {
  journal::Method<journal::ExternalVesselGetPosition> m{
      {plugin,
       vessel_guid,
       time},
      {position}};
  if (plugin == nullptr) {
    return m.Return(
        MakeStatus(Error::INVALID_ARGUMENT, "|plugin| must not be null"));
  }
  Instant const t = FromGameTime(*plugin, time);
  if (!plugin->HasVessel(vessel_guid)) {
    return m.Return(MakeStatus(
        Error::NOT_FOUND,
        absl::StrCat("No vessel with GUID ", vessel_guid)));
  }
  auto const& vessel = *plugin->GetVessel(vessel_guid);
  auto const& trajectory = vessel.psychohistory();
  if (t < trajectory.t_min() || t > trajectory.t_max()) {
    return m.Return(MakeStatus(
        Error::OUT_OF_RANGE,
        (std::stringstream{}
         << "|time| " << t << " does not lie within the domain ["
         << trajectory.t_min() << ", " << trajectory.t_max()
         << "] of the psychohistory of " << vessel.ShortDebugString()).str()));
  }
  auto const from_solar_system_barycentre =
      plugin->renderer().BarycentricToWorldSun(plugin->PlanetariumRotation())(
          trajectory.EvaluatePosition(t) - Barycentric::origin);
  *position = ToXYZ(from_solar_system_barycentre.coordinates() / Metre);
  return m.Return(OK());
}

}  // namespace interface
}  // namespace principia
