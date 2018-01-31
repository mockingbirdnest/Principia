
#include "ksp_plugin/interface.hpp"

#include "base/array.hpp"
#include "base/status.hpp"
#include "base/status_or.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"

namespace principia {
namespace interface {

using base::Error;
using base::Status;
using base::StatusOr;
using base::UniqueBytes;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::Vessel;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DiscreteTrajectory;

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

int principia__ExternalFlowBodyCentred(
    Plugin const* const plugin,
    int const central_body_index,
    QP const world_body_centred_initial_degrees_of_freedom,
    double const t_initial,
    double const t_final,
    QP* const world_body_centred_final_degrees_of_freedom,
    char const** const error_message) {
  journal::Method<journal::ExternalFlowBodyCentred> m{
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
  return error_code(Error::UNIMPLEMENTED);
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
               " out of range, vessel " + vessel.ShortDebugString() + " has " +
               flight_plan.number_of_manœuvres() + u8" planned manœuvres"),
        error_message));
  }
  // The index of the coast segment following the desired manœuvre.
  int const segment_index = manoeuvre_index * 2 + 3;
  if (segment_index >= flight_plan.number_of_segments()) {
    return m.Return(set_error(
        Status(Error::FAILED_PRECONDITION,
               u8"A singularity occurs within manœuvre " + manoeuvre_index +
                   " of " + vessel.ShortDebugString()),
        error_message))
  }
  DiscreteTrajectory<Barycentric>::Iterator coast_begin;
  DiscreteTrajectory<Barycentric>::Iterator coast_end;
  flight_plan.GetSegment(segment_index, begin, end);
  auto const body_centred_inertial =
      plugin->NewBodyCentredNonRotatingNavigationFrame(central_body_index);
  DiscreteTrajectory<Navigation> body_centred_inertial_coast;
  for (auto it = coast_begin; it != coast_end; ++it) {
    body_centred_inertial_coast.Append(
        it.time(),
        body_centred_inertial->ToThisFrameAtTime(it.time())(
            it.degrees_of_freedom()));
  }
  plugin->renderer().WorldToBarycentric(plugin->PlanetariumRotation())(
      FromXYZ<Displacement<World>>(world_body_centred_reference_position));
  DiscreteTrajectory<Navigation> immobile_reference_point;
  immobile_reference_point.Append(body_centred_inertial_coast.Begin().time(), )
}

}  // namespace interface
}  // namespace principia
