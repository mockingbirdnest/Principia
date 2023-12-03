#include "ksp_plugin/interface.hpp"

#include "base/push_pull_callback.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/apsides.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::base::_push_pull_callback;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_apsides;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

QP __cdecl principia__CollisionDeleteExecutor(
    PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>** const
        executor) {
  journal::Method<journal::CollisionDeleteExecutor> m{{executor}, {executor}};
  CHECK_NOTNULL(executor);
  auto const collision_degrees_of_freedom = (*executor)->get();
  {
    TakeOwnership(executor);
  }
  return ToQP(collision_degrees_of_freedom);
}

bool __cdecl principia__CollisionGetLatitudeLongitude(
    PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>* const
        executor,
    double* const latitude_in_degrees,
    double* const longitude_in_degrees) {
  journal::Method<journal::CollisionGetLatitudeLongitude> m{
      {executor},
      {latitude_in_degrees,
       longitude_in_degrees}};
  CHECK_NOTNULL(executor);

  Angle latitude;
  Angle longitude;
  bool const more = executor->callback().Pull(latitude, longitude);
  *latitude_in_degrees = latitude / Degree;
  *longitude_in_degrees = longitude / Degree;

  return m.Return(more);
}

PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>*
__cdecl principia__CollisionNewExecutor(
    Plugin const* const plugin,
    int const celestial_index,
    XYZ const sun_world_position,
    char const* const vessel_guid,
    double const apoapside_time,
    double const periapside_time) {
  journal::Method<journal::CollisionNewExecutor> m{
      {plugin,
       celestial_index,
       sun_world_position,
       vessel_guid,
       apoapside_time,
       periapside_time}};
  CHECK_NOTNULL(plugin);
  not_null<Vessel*> const vessel = plugin->GetVessel(vessel_guid);
  auto const vessel_prediction = vessel->prediction();
  auto const begin =
      vessel_prediction->lower_bound(FromGameTime(*plugin, apoapside_time));
  auto const end =
      vessel_prediction->lower_bound(FromGameTime(*plugin, periapside_time));

  auto task = [begin,
               celestial_index,
               end,
               plugin,
               sun_world_position =
                   FromXYZ<Position<World>>(sun_world_position),
               vessel_prediction](
                  std::function<Length(Angle const& latitude,
                                       Angle const& longitude)> const& radius) {
    return plugin->ComputeAndRenderCollision(celestial_index,
                                             *vessel_prediction,
                                             begin,
                                             end,
                                             sun_world_position,
                                             radius);
  };

  return m.Return(
      new PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>(
          std::move(task)));
}

void __cdecl principia__CollisionSetRadius(
    PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>* const
        executor,
    double const radius) {
  journal::Method<journal::CollisionSetRadius> m{
      {executor,
       radius}};
  CHECK_NOTNULL(executor);
  executor->callback().Push(radius * Metre);
}

}  // namespace interface
}  // namespace principia
