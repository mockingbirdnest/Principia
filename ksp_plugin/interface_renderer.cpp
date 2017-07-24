
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

// Returns the frame last set by |plugin->SetPlottingFrame|.  No transfer of
// ownership.  The returned pointer is never null.
NavigationFrame const* principia__GetPlottingFrame(Plugin const* const plugin) {
  journal::Method<journal::GetPlottingFrame> m({plugin});
  return m.Return(GetRenderer(plugin).GetPlottingFrame());
}

Iterator* principia__RenderedPrediction(Plugin* const plugin,
                                        char const* const vessel_guid,
                                        XYZ const sun_world_position) {
  journal::Method<journal::RenderedPrediction> m({plugin,
                                                  vessel_guid,
                                                  sun_world_position});
  CHECK_NOTNULL(plugin);
  auto const& prediction = plugin->GetVessel(vessel_guid)->prediction();
  auto rendered_trajectory =
      GetRenderer(plugin).RenderBarycentricTrajectoryInWorld(
          plugin->CurrentTime(),
          prediction.Begin(),
          prediction.End(),
          FromXYZ<Position<World>>(sun_world_position),
          plugin->PlanetariumRotation());
  return m.Return(new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_trajectory),
      plugin));
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
  auto const& prediction = plugin->GetVessel(vessel_guid)->prediction();
  std::unique_ptr<DiscreteTrajectory<World>> rendered_apoapsides;
  std::unique_ptr<DiscreteTrajectory<World>> rendered_periapsides;
  plugin->ComputeAndRenderApsides(celestial_index,
                                  prediction.Begin(),
                                  prediction.End(),
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
  auto const& prediction = plugin->GetVessel(vessel_guid)->prediction();
  std::unique_ptr<DiscreteTrajectory<World>> rendered_closest_approaches;
  plugin->ComputeAndRenderClosestApproaches(
      prediction.Begin(),
      prediction.End(),
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
  auto const& prediction = plugin->GetVessel(vessel_guid)->prediction();
  std::unique_ptr<DiscreteTrajectory<World>> rendered_ascending;
  std::unique_ptr<DiscreteTrajectory<World>> rendered_descending;
  plugin->ComputeAndRenderNodes(prediction.Begin(),
                                prediction.End(),
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

Iterator* principia__RenderedVesselTrajectory(Plugin const* const plugin,
                                              char const* const vessel_guid,
                                              XYZ const sun_world_position) {
  journal::Method<journal::RenderedVesselTrajectory> m({plugin,
                                                        vessel_guid,
                                                        sun_world_position});
  CHECK_NOTNULL(plugin);
  auto const& psychohistory = plugin->GetVessel(vessel_guid)->psychohistory();
  auto rendered_trajectory =
      GetRenderer(plugin).RenderBarycentricTrajectoryInWorld(
          plugin->CurrentTime(),
          psychohistory.Begin(),
          psychohistory.End(),
          FromXYZ<Position<World>>(sun_world_position),
          plugin->PlanetariumRotation());
  return m.Return(new TypedIterator<DiscreteTrajectory<World>>(
      std::move(rendered_trajectory),
      plugin));
}

// |navigation_frame| must not be null.  No transfer of ownership of
// |*navigation_frame|, takes ownership of |**navigation_frame|, nulls
// |*navigation_frame|.
void principia__SetPlottingFrame(Plugin* const plugin,
                                 NavigationFrame** const navigation_frame) {
  journal::Method<journal::SetPlottingFrame> m({plugin, navigation_frame},
                                               {navigation_frame});
  GetRenderer(plugin).SetPlottingFrame(TakeOwnership(navigation_frame));
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
