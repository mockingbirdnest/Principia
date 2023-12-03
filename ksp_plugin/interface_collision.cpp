#include "ksp_plugin/interface.hpp"

#include "base/push_pull_callback.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.

namespace principia {
namespace interface {

using namespace principia::base::_push_pull_callback;
using namespace principia::journal::_method;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_apsides;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_quantities;

PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>*
__cdecl principia__CollisionNewExecutor(
    Plugin const* const plugin,
    int const celestial_index,
    char const* const vessel_guid,
    double const apoapside_time,
    double const periapside_time) {
  journal::Method<journal::CollisionNewExecutor> m{
      {plugin,
       celestial_index,
       vessel_guid,
       apoapside_time,
       periapside_time}};
  CHECK_NOTNULL(plugin);
  not_null<Vessel*> const vessel = plugin->GetVessel(vessel_guid);
  auto const vessel_prediction = vessel->prediction();
  //TODO(phl)lower_bound
  auto const begin =
      vessel_prediction->find(FromGameTime(*plugin, apoapside_time));
  auto const end =
      vessel_prediction->find(FromGameTime(*plugin, periapside_time));

  auto task = [begin, celestial_index, end, plugin, vessel_prediction](
                  std::function<Length(Angle const& latitude,
                                       Angle const& longitude)> const& radius) {
    return plugin->ComputeAndRenderCollision(
        celestial_index,
        *vessel_prediction,
        begin,
        end,
        /*sun_world_position=*/Position<World>(),//TODO(phl)nonono
        radius);
  };

  return m.Return(
      new PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>(
          std::move(task)));
}

void __cdecl principia__CollisionDeleteExecutor(
    PushPullExecutor<DegreesOfFreedom<World>, Length, Angle, Angle>
        const** const executor) {}

bool __cdecl principia__CollisionGetLatitudeLongitude(double* const latitude,
                                                      double* const longitude) {
}

void __cdecl principia__CollisionSetRadius(double const radius) {}

}  // namespace interface
}  // namespace principia
