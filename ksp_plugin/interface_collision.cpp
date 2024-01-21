#include "ksp_plugin/interface.hpp"

#include <memory>
#include <utility>

#include "base/push_pull_callback.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::base::_push_pull_callback;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

namespace {

template<typename TrajectoryLike>
not_null<std::unique_ptr<
    PushPullExecutor<std::optional<DiscreteTrajectory<World>::value_type>,
                     Length, Angle, Angle>>>
NewExecutor(Plugin const* const plugin,
            int const celestial_index,
            XYZ const sun_world_position,
            int const max_points,
            TrajectoryLike const& vessel_trajectory) {
  CHECK_NOTNULL(plugin);

  auto task = [celestial_index,
               max_points,
               plugin,
               sun_world_position =
                   FromXYZ<Position<World>>(sun_world_position),
               &vessel_trajectory](
                  std::function<Length(Angle const& latitude,
                                       Angle const& longitude)> const& radius) {
    return plugin->ComputeAndRenderFirstCollision(celestial_index,
                                                  vessel_trajectory,
                                                  vessel_trajectory.begin(),
                                                  vessel_trajectory.end(),
                                                  sun_world_position,
                                                  max_points,
                                                  radius);
  };

  return make_not_null_unique<
      PushPullExecutor<std::optional<DiscreteTrajectory<World>::value_type>,
                       Length, Angle, Angle>>(std::move(task));
}


}  // namespace

bool __cdecl principia__CollisionDeleteExecutor(
    Plugin const* const plugin,
    PushPullExecutor<std::optional<DiscreteTrajectory<World>::value_type>,
                     Length, Angle, Angle>** const executor,
    TQP* const collision) {
  journal::Method<journal::CollisionDeleteExecutor> m{{plugin, executor},
                                                      {executor, collision}};
  CHECK_NOTNULL(executor);
  auto const maybe_collision = (*executor)->get();
  {
    TakeOwnership(executor);
  }
  if (maybe_collision.has_value()) {
    *collision = TQP{.t = ToGameTime(*plugin, maybe_collision->time),
                     .qp = ToQP(maybe_collision->degrees_of_freedom)};
    return m.Return(true);
  } else {
    return m.Return(false);
  }
}

bool __cdecl principia__CollisionGetLatitudeLongitude(
    PushPullExecutor<std::optional<DiscreteTrajectory<World>::value_type>,
                     Length, Angle, Angle>* const executor,
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

PushPullExecutor<
    std::optional<DiscreteTrajectory<World>::value_type>,
    Length, Angle, Angle>* __cdecl principia__CollisionNewFlightPlanExecutor(
    Plugin const* const plugin,
    int const celestial_index,
    XYZ const sun_world_position,
    int const max_points,
    char const* const vessel_guid) {
  journal::Method<journal::CollisionNewFlightPlanExecutor> m{
      {plugin,
       celestial_index,
       sun_world_position,
       max_points,
       vessel_guid}};
  CHECK_NOTNULL(plugin);
  auto& flight_plan = GetFlightPlan(*plugin, vessel_guid);
  return m.Return(NewExecutor(plugin,
                              celestial_index,
                              sun_world_position,
                              max_points,
                              flight_plan.GetAllSegments())
                      .release());
}

PushPullExecutor<
    std::optional<DiscreteTrajectory<World>::value_type>,
    Length, Angle, Angle>* __cdecl principia__CollisionNewPredictionExecutor(
    Plugin const* const plugin,
    int const celestial_index,
    XYZ const sun_world_position,
    int const max_points,
    char const* const vessel_guid) {
  journal::Method<journal::CollisionNewPredictionExecutor> m{
      {plugin,
       celestial_index,
       sun_world_position,
       max_points,
       vessel_guid}};
  CHECK_NOTNULL(plugin);
  not_null<Vessel*> const vessel = plugin->GetVessel(vessel_guid);
  return m.Return(NewExecutor(plugin,
                              celestial_index,
                              sun_world_position,
                              max_points,
                              *vessel->prediction())
                      .release());
}

void __cdecl principia__CollisionSetRadius(
    PushPullExecutor<std::optional<DiscreteTrajectory<World>::value_type>,
                     Length, Angle, Angle>* const executor,
    double const radius) {
  journal::Method<journal::CollisionSetRadius> m{
      {executor,
       radius}};
  CHECK_NOTNULL(executor);
  executor->callback().Push(radius * Metre);
  return m.Return();
}

}  // namespace interface
}  // namespace principia
