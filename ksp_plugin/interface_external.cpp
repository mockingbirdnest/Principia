
#include "ksp_plugin/interface.hpp"

#include "base/array.hpp"
#include "base/status.hpp"
#include "base/status_or.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "physics/apsides.hpp"

namespace principia {
namespace interface {

using base::Error;
using base::Status;
using base::StatusOr;
using base::UniqueBytes;
using geometry::AngularVelocity;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::Vessel;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::ComputeApsides;
using physics::DiscreteTrajectory;
using physics::RigidMotion;
using physics::RigidTransformation;

namespace {

char const* release_string(const std::string& s) {
  UniqueBytes allocated_string(s.size() + 1);
  std::memcpy(allocated_string.data.get(), s.data(), s.size() + 1);
  return reinterpret_cast<char const*>(allocated_string.data.release());
}

int set_error(Status status, char const** error_message) {
  *error_message = release_string(status.message());
  return static_cast<int>(status.error());
}

}  // namespace

int principia__ExternalFlowFreefall(
    Plugin const* const plugin,
    int const central_body_index,
    QP const world_body_centred_initial_degrees_of_freedom,
    double const t_initial,
    double const t_final,
    QP* const world_body_centred_final_degrees_of_freedom,
    char const** const error_message) {
  journal::Method<journal::ExternalFlowFreefall> m{
      {plugin,
       central_body_index,
       world_body_centred_initial_degrees_of_freedom,
       t_initial,
       t_final},
      {world_body_centred_final_degrees_of_freedom, error_message}};
  if (plugin == nullptr) {
    return m.Return(
        set_error(Status(Error::INVALID_ARGUMENT, "|plugin| must not be null"),
                  error_message));
  }
  return set_error(Status(Error::UNIMPLEMENTED,
                          "|ExternalFlowFreefall| is not yet implemented"),
                   error_message);
}

int principia__ExternalGetNearestPlannedCoastDegreesOfFreedom(
    Plugin const* const plugin,
    int const central_body_index,
    char const* const vessel_guid,
    int manoeuvre_index,
    XYZ world_body_centred_reference_position,
    QP* world_body_centred_nearest_degrees_of_freedom,
    char const** const error_message) {
  journal::Method<journal::ExternalGetNearestPlannedCoastDegreesOfFreedom> m{
      {plugin,
       central_body_index,
       vessel_guid,
       manoeuvre_index,
       world_body_centred_reference_position},
      {world_body_centred_nearest_degrees_of_freedom, error_message}};
  if (plugin == nullptr) {
    return m.Return(
        set_error(Status(Error::INVALID_ARGUMENT, "|plugin| must not be null"),
                  error_message));
  }
  if (manoeuvre_index < 0) {
    return m.Return(set_error(Status(Error::INVALID_ARGUMENT,
                                     "Invalid negative |manoeuvre_index|" +
                                         std::to_string(manoeuvre_index)),
                              error_message));
  }
  if (!plugin->HasCelestial(central_body_index)) {
    return m.Return(set_error(
        Status(Error::NOT_FOUND,
               "No celestial with index " + std::to_string(central_body_index)),
        error_message));
  }
  if (!plugin->HasVessel(vessel_guid)) {
    return m.Return(
        set_error(Status(Error::NOT_FOUND,
                         "No vessel with GUID " + std::string(vessel_guid)),
                  error_message));
  }
  Vessel const& vessel = *plugin->GetVessel(vessel_guid);
  if (!vessel.has_flight_plan()) {
    return m.Return(set_error(
        Status(Error::FAILED_PRECONDITION,
               "Vessel " + vessel.ShortDebugString() + " has no flight plan"),
        error_message));
  }
  FlightPlan const& flight_plan = vessel.flight_plan();
  if (manoeuvre_index >= flight_plan.number_of_manœuvres()) {
    return m.Return(set_error(
        Status(Error::OUT_OF_RANGE,
               "|manoeuvre_index| " + std::to_string(manoeuvre_index) +
                   " out of range, vessel " + vessel.ShortDebugString() +
                   " has " + std::to_string(flight_plan.number_of_manœuvres()) +
                   u8" planned manœuvres"),
        error_message));
  }
  // The index of the coast segment following the desired manœuvre.
  int const segment_index = manoeuvre_index * 2 + 3;
  if (segment_index >= flight_plan.number_of_segments()) {
    return m.Return(set_error(Status(Error::FAILED_PRECONDITION,
                                     u8"A singularity occurs within manœuvre " +
                                         std::to_string(manoeuvre_index) +
                                         " of " + vessel.ShortDebugString()),
                              error_message));
  }
  DiscreteTrajectory<Barycentric>::Iterator coast_begin;
  DiscreteTrajectory<Barycentric>::Iterator coast_end;
  flight_plan.GetSegment(segment_index, coast_begin, coast_end);
  auto const body_centred_inertial =
      plugin->NewBodyCentredNonRotatingNavigationFrame(central_body_index);
  DiscreteTrajectory<Navigation> coast;
  for (auto it = coast_begin; it != coast_end; ++it) {
    coast.Append(it.time(),
                 body_centred_inertial->ToThisFrameAtTime(it.time())(
                     it.degrees_of_freedom()));
  }

  Instant const current_time = plugin->CurrentTime();
  // The given |World| position and requested |World| degrees of freedom are
  // body-centred inertial, so |body_centred_inertial| up to an orthogonal map
  // to world coordinates.  Do the conversion directly.
  // NOTE(eggrobin): it is correct to use the orthogonal map at |current_time|,
  // because |body_centred_inertial| does not rotate with respect to
  // |Barycentric|, so the orthogonal map does not depend on time.
  RigidMotion<Navigation, World> to_world_body_centred_inertial(
      RigidTransformation<Navigation, World>(
          Navigation::origin,
          World::origin,
          plugin->renderer().BarycentricToWorld(plugin->PlanetariumRotation()) *
              body_centred_inertial->FromThisFrameAtTime(
                  current_time).orthogonal_map()),
      AngularVelocity<Navigation>{},
      Velocity<Navigation>{});
  auto const from_world_body_centred_inertial =
      to_world_body_centred_inertial.Inverse();
  Position<Navigation> reference_position =
      from_world_body_centred_inertial.rigid_transformation()(
          FromXYZ<Position<World>>(world_body_centred_reference_position));
  DiscreteTrajectory<Navigation> immobile_reference;
  immobile_reference.Append(coast.Begin().time(),
                            {reference_position, Velocity<Navigation>{}});
  if (coast.Begin() !=
      coast.last()) {
    immobile_reference.Append(coast.last().time(),
                              {reference_position, Velocity<Navigation>{}});
  }
  DiscreteTrajectory<Navigation> apoapsides;
  DiscreteTrajectory<Navigation> periapsides;
  ComputeApsides(/*reference=*/immobile_reference,
                 coast.Begin(),
                 coast.End(),
                 apoapsides,
                 periapsides);
  if (periapsides.Empty()) {
    bool const coasting_away =
        (coast.Begin().degrees_of_freedom().position() -
         reference_position).Norm²() <
        (coast.last().degrees_of_freedom().position() -
         reference_position).Norm²();
    *world_body_centred_nearest_degrees_of_freedom =
        ToQP(to_world_body_centred_inertial(
            coasting_away ? coast.Begin().degrees_of_freedom()
                          : coast.last().degrees_of_freedom()));
  } else {
    *world_body_centred_nearest_degrees_of_freedom =
        ToQP(to_world_body_centred_inertial(
            periapsides.Begin().degrees_of_freedom()));
  }
  return m.Return(set_error(Status::OK, error_message));
}

}  // namespace interface
}  // namespace principia
