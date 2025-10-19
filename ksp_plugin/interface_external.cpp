#include "ksp_plugin/interface.hpp"

#include <limits>
#include <string>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "base/array.hpp"
#include "geometry/frame.hpp"
#include "geometry/r3_element.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace interface {

using namespace principia::base::_array;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_r3_element;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_apsides;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_rigid_motion;

namespace {

Status* Unimplemented() {
  static Status* const unimplemented =
      ToNewStatus(absl::UnimplementedError("API entry point was deprecated"));
  return unimplemented;
}

}  // namespace

Status* __cdecl principia__ExternalCelestialGetPosition(
    Plugin const* const plugin,
    int const body_index,
    double const time,
    XYZ* const position) {
  journal::Method<journal::ExternalCelestialGetPosition> m{
      {plugin,
       body_index,
       time},
      {position}};
  return m.Return(Unimplemented());
}

Status* __cdecl principia__ExternalCelestialGetSurfacePosition(
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
  return m.Return(Unimplemented());
}

Status* __cdecl principia__ExternalFlowFreefall(
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
  return m.Return(Unimplemented());
}

Status* __cdecl principia__ExternalGeopotentialGetCoefficient(
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
  return m.Return(Unimplemented());
}

Status* __cdecl principia__ExternalGeopotentialGetReferenceRadius(
    Plugin const* const plugin,
    int const body_index,
    double* const reference_radius) {
  journal::Method<journal::ExternalGeopotentialGetReferenceRadius> m{
      {plugin,
       body_index},
      {reference_radius}};
  return m.Return(Unimplemented());
}

Status* __cdecl principia__ExternalGetNearestPlannedCoastDegreesOfFreedom(
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
  return m.Return(Unimplemented());
}

Status* __cdecl principia__ExternalVesselGetPosition(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const time,
    XYZ* const position) {
  journal::Method<journal::ExternalVesselGetPosition> m{
      {plugin,
       vessel_guid,
       time},
      {position}};
  return m.Return(Unimplemented());
}

}  // namespace interface
}  // namespace principia
