#include "ksp_plugin/renderer.hpp"

#include "geometry/grassmann.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using base::make_not_null_unique;
using geometry::AngularVelocity;
using geometry::Vector;
using geometry::Velocity;
using physics::BodyCentredBodyDirectionDynamicFrame;
using physics::ComputeApsides;
using physics::ComputeNodes;
using physics::DegreesOfFreedom;
using physics::RigidTransformation;

Renderer::Renderer(not_null<Celestial const*> const sun,
                   not_null<std::unique_ptr<NavigationFrame>> plotting_frame)
    : sun_(sun),
      plotting_frame_(std::move(plotting_frame)) {}

void Renderer::SetPlottingFrame(
    not_null<std::unique_ptr<NavigationFrame>> plotting_frame) {
  plotting_frame_ = std::move(plotting_frame);
}

not_null<NavigationFrame const*> Renderer::GetPlottingFrame() const {
  return target_ ? target_->target_frame.get()
                 : plotting_frame_.get();
}

void Renderer::SetTargetVessel(
    not_null<Vessel*> const vessel,
    not_null<Celestial const*> const celestial,
    not_null<Ephemeris<Barycentric> const*> const ephemeris) {
  CHECK(!vessel->prediction().Empty());
  if (!target_ ||
      target_->vessel != vessel ||
      target_->celestial != celestial) {
    target_.emplace(vessel, celestial, ephemeris);
  }
}

void Renderer::ClearTargetVessel() {
  target_ = std::experimental::nullopt;
}

bool Renderer::HasTargetVessel() const {
  return (bool)target_;
}

Vessel const & Renderer::GetTargetVessel() const {
  CHECK(target_);
  return *target_->vessel;
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Renderer::RenderBarycentricTrajectoryInWorld(
    Instant const& time,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  auto const trajectory_in_plotting_frame =
      RenderBarycentricTrajectoryInPlotting(begin, end);
  auto trajectory_in_world =
      RenderPlottingTrajectoryInWorld(time,
                                      trajectory_in_plotting_frame->Begin(),
                                      trajectory_in_plotting_frame->End(),
                                      sun_world_position,
                                      planetarium_rotation);
  return trajectory_in_world;
}

not_null<std::unique_ptr<DiscreteTrajectory<Navigation>>>
Renderer::RenderBarycentricTrajectoryInPlotting(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end) const {
  CHECK(!target_ || begin == end ||
        (!target_->vessel->prediction().Empty() &&
         end.time() <= target_->vessel->prediction().last().time()));

  auto trajectory = make_not_null_unique<DiscreteTrajectory<Navigation>>();
  for (auto it = begin; it != end; ++it) {
    if (target_) {
      if (it.time() < target_->vessel->prediction().t_min()) {
        continue;
      } else if (it.time() > target_->vessel->prediction().t_max()) {
        break;
      }
    }
    trajectory->Append(
        it.time(),
        BarycentricToPlotting(it.time())(it.degrees_of_freedom()));
  }
  VLOG(1) << "Returning a " << trajectory->Size() << "-point trajectory";
  return trajectory;
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Renderer::RenderPlottingTrajectoryInWorld(
    Instant const& time,
    DiscreteTrajectory<Navigation>::Iterator const& begin,
    DiscreteTrajectory<Navigation>::Iterator const& end,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<World>>();

  RigidMotion<Navigation, World> const
      from_plotting_frame_to_world_at_current_time =
          PlottingToWorld(time, sun_world_position, planetarium_rotation);
  for (auto it = begin; it != end; ++it) {
    DegreesOfFreedom<Navigation> const& navigation_degrees_of_freedom =
        it.degrees_of_freedom();
    DegreesOfFreedom<World> const world_degrees_of_freedom =
        from_plotting_frame_to_world_at_current_time(
            navigation_degrees_of_freedom);
    trajectory->Append(it.time(), world_degrees_of_freedom);
  }
  VLOG(1) << "Returning a " << trajectory->Size() << "-point trajectory";
  return trajectory;
}

RigidMotion<Barycentric, Navigation> Renderer::BarycentricToPlotting(
    Instant const& time) const {
  return GetPlottingFrame()->ToThisFrameAtTime(time);
}

RigidMotion<Barycentric, World> Renderer::BarycentricToWorld(
    Instant const& time,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return RigidMotion<Barycentric, World>(
      RigidTransformation<Barycentric, World>(
          sun_->current_position(time),
          sun_world_position,
          BarycentricToWorld(planetarium_rotation)),
      AngularVelocity<Barycentric>{},
      Velocity<Barycentric>{});
}

OrthogonalMap<Barycentric, World> Renderer::BarycentricToWorld(
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return OrthogonalMap<WorldSun, World>::Identity() *
         BarycentricToWorldSun(planetarium_rotation);
}

OrthogonalMap<Barycentric, WorldSun> Renderer::BarycentricToWorldSun(
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return sun_looking_glass.Inverse().Forget() * planetarium_rotation.Forget();
}

OrthogonalMap<Frenet<Navigation>, World> Renderer::FrenetToWorld(
    Vessel const& vessel,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  auto const last = vessel.psychohistory().last();
  Instant const& time = last.time();
  DegreesOfFreedom<Barycentric> const& barycentric_degrees_of_freedom =
      last.degrees_of_freedom();
  DegreesOfFreedom<Navigation> const plotting_frame_degrees_of_freedom =
      BarycentricToPlotting(time)(barycentric_degrees_of_freedom);
  Rotation<Frenet<Navigation>, Navigation> const
      from_frenet_frame_to_plotting_frame =
          GetPlottingFrame()->FrenetFrame(time,
                                          plotting_frame_degrees_of_freedom);

  return PlottingToWorld(time, planetarium_rotation) *
         from_frenet_frame_to_plotting_frame.Forget();
}

OrthogonalMap<Navigation, Barycentric> Renderer::PlottingToBarycentric(
    Instant const& time) const {
  return GetPlottingFrame()->FromThisFrameAtTime(time).orthogonal_map();
}

RigidMotion<Navigation, World> Renderer::PlottingToWorld(
    Instant const& time,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return BarycentricToWorld(time, sun_world_position, planetarium_rotation) *
         GetPlottingFrame()->FromThisFrameAtTime(time);
}

OrthogonalMap<Navigation, World> Renderer::PlottingToWorld(
    Instant const& time,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return BarycentricToWorld(planetarium_rotation) *
         PlottingToBarycentric(time);
}

RigidMotion<World, Barycentric> Renderer::WorldToBarycentric(
    Instant const& time,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return BarycentricToWorld(time, sun_world_position, planetarium_rotation)
      .Inverse();
}

OrthogonalMap<World, Barycentric> Renderer::WorldToBarycentric(
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return BarycentricToWorld(planetarium_rotation).Inverse();
}

Renderer::Target::Target(
    not_null<Vessel*> const vessel,
    not_null<Celestial const*> const celestial,
    not_null<Ephemeris<Barycentric> const*> const ephemeris)
    : vessel(vessel),
      celestial(celestial),
      target_frame(
          make_not_null_unique<
              BodyCentredBodyDirectionDynamicFrame<Barycentric, Navigation>>(
              ephemeris,
              [this]() -> auto& { return this->vessel->prediction(); },
              celestial->body())) {}

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
