
#include "ksp_plugin/interface.hpp"

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/renderer.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace interface {

using ksp_plugin::Renderer;
using ksp_plugin::TypedIterator;
using physics::DiscreteTrajectory;

namespace {

Renderer& GetRenderer(Plugin* const plugin) {
  return CHECK_NOTNULL(plugin)->renderer();
}

Renderer const& GetRenderer(Plugin const* const plugin) {
  return CHECK_NOTNULL(plugin)->renderer();
}

}  // namespace

void principia__ClearTargetVessel(Plugin* const plugin) {
  journal::Method<journal::ClearTargetVessel> m({plugin});
  GetRenderer(plugin).ClearTargetVessel();
  return m.Return();
}

void principia__RenderedPredictionApsides(Plugin const* const plugin,
                                          char const* const vessel_guid,
                                          int const celestial_index,
                                          XYZ const sun_world_position,
                                          Iterator** const apoapsides,
                                          Iterator** const periapsides) {
  journal::Method<journal::RenderedPredictionApsides> m(
      {plugin, vessel_guid, celestial_index, sun_world_position},
      {apoapsides, periapsides});
  CHECK_NOTNULL(plugin);
  auto const prediction = plugin->GetVessel(vessel_guid)->prediction();
  std::unique_ptr<DiscreteTrajectory<World>> rendered_apoapsides;
  std::unique_ptr<DiscreteTrajectory<World>> rendered_periapsides;
  plugin->ComputeAndRenderApsides(celestial_index,
                                  prediction->Fork(),
                                  prediction->End(),
                                  FromXYZ<Position<World>>(sun_world_position),
                                  rendered_apoapsides,
                                  rendered_periapsides);
  *apoapsides = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_apoapsides)),
      plugin);
  *periapsides = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_periapsides)),
      plugin);
  return m.Return();
}

void principia__RenderedPredictionClosestApproaches(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    Iterator** const closest_approaches) {
  journal::Method<journal::RenderedPredictionClosestApproaches> m(
      {plugin, vessel_guid, sun_world_position},
      {closest_approaches});
  CHECK_NOTNULL(plugin);
  auto const prediction = plugin->GetVessel(vessel_guid)->prediction();
  std::unique_ptr<DiscreteTrajectory<World>> rendered_closest_approaches;
  plugin->ComputeAndRenderClosestApproaches(
      prediction->Fork(),
      prediction->End(),
      FromXYZ<Position<World>>(sun_world_position),
      rendered_closest_approaches);
  *closest_approaches = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_closest_approaches)),
      plugin);
  return m.Return();
}

void principia__RenderedPredictionNodes(Plugin const* const plugin,
                                        char const* const vessel_guid,
                                        XYZ const sun_world_position,
                                        Iterator** const ascending,
                                        Iterator** const descending) {
  journal::Method<journal::RenderedPredictionNodes> m(
      {plugin, vessel_guid, sun_world_position},
      {ascending, descending});
  CHECK_NOTNULL(plugin);
  auto const prediction = plugin->GetVessel(vessel_guid)->prediction();
  std::unique_ptr<DiscreteTrajectory<World>> rendered_ascending;
  std::unique_ptr<DiscreteTrajectory<World>> rendered_descending;
  plugin->ComputeAndRenderNodes(prediction->Fork(),
                                prediction->End(),
                                FromXYZ<Position<World>>(sun_world_position),
                                rendered_ascending,
                                rendered_descending);
  *ascending = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_ascending)),
      plugin);
  *descending = new TypedIterator<DiscreteTrajectory<World>>(
      check_not_null(std::move(rendered_descending)),
      plugin);
  return m.Return();
}

// Calls |plugin| to create a |NavigationFrame| using the given |parameters|,
// sets it as the current plotting frame.
void principia__SetPlottingFrame(Plugin* const plugin,
                                 NavigationFrameParameters const parameters) {
  journal::Method<journal::SetPlottingFrame> m({plugin, parameters});
  auto navigation_frame = NewNavigationFrame(*plugin, parameters);
  GetRenderer(plugin).SetPlottingFrame(std::move(navigation_frame));
  return m.Return();
}

void principia__SetTargetVessel(Plugin* const plugin,
                                char const* const vessel_guid,
                                int const reference_body_index) {
  journal::Method<journal::SetTargetVessel> m(
      {plugin, vessel_guid, reference_body_index});
  CHECK_NOTNULL(plugin);
  plugin->SetTargetVessel(vessel_guid, reference_body_index);
  return m.Return();
}

}  // namespace interface
}  // namespace principia
