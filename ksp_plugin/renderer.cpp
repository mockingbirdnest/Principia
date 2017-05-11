#include "ksp_plugin/renderer.hpp"

#include "geometry/grassmann.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_body_direction_dynamic_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/rigid_motion.hpp"

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
using physics::RigidMotion;

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
    Instant const& time,
    not_null<Vessel*> const vessel,
    not_null<Celestial const*> const celestial,
    not_null<Ephemeris<Barycentric> const*> const ephemeris) {
  if (!target_ || target_->vessel != vessel ||
      target_->celestial != celestial) {
    target_.emplace(vessel, celestial, ephemeris);
  }
  // Make sure that the current time is covered by the prediction.
  if (time > target_->vessel->prediction().t_max()) {
    target_->vessel->UpdatePrediction(time + prediction_length_);
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
  auto const trajectory_in_navigation =
      RenderBarycentricTrajectoryInNavigation(time, begin, end);
  auto trajectory_in_world =
      RenderNavigationTrajectoryInWorld(time,
                                        trajectory_in_navigation->Begin(),
                                        trajectory_in_navigation->End(),
                                        sun_world_position,
                                        planetarium_rotation);
  return trajectory_in_world;
}

not_null<std::unique_ptr<DiscreteTrajectory<Navigation>>>
Renderer::RenderBarycentricTrajectoryInNavigation(
    Instant const& time,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end) const {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<Navigation>>();

  NavigationFrame const& plotting_frame = *GetPlottingFrame();

  if (target_ && !begin.trajectory()->Empty() &&
      (target_->vessel->prediction().Empty() ||
       begin.trajectory()->last().time() >
           target_->vessel->prediction().last().time())) {
    // NOTE(egg): this is an ugly hack to try to get a long enough trajectory
    // while retaining a timeout.
    auto parameters = target_->vessel->prediction_adaptive_step_parameters();
    parameters.set_max_steps(begin.trajectory()->Size());
    target_->vessel->set_prediction_adaptive_step_parameters(parameters);
    target_->vessel->UpdatePrediction(time + prediction_length_);
  }

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
        plotting_frame.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
  }
  VLOG(1) << "Returning a " << trajectory->Size() << "-point trajectory";
  return trajectory;
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Renderer::RenderNavigationTrajectoryInWorld(
    Instant const& time,
    DiscreteTrajectory<Navigation>::Iterator const& begin,
    DiscreteTrajectory<Navigation>::Iterator const& end,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<World>>();

  NavigationFrame const& plotting_frame = *GetPlottingFrame();

  RigidMotion<Navigation, World> from_navigation_frame_to_world_at_current_time(
      /*rigid_transformation=*/
          BarycentricToWorld(time, sun_world_position, planetarium_rotation) *
          plotting_frame.FromThisFrameAtTime(time).rigid_transformation(),
      AngularVelocity<Navigation>{},
      Velocity<Navigation>{});
  for (auto it = begin; it != end; ++it) {
    DegreesOfFreedom<Navigation> const& navigation_degrees_of_freedom =
        it.degrees_of_freedom();
    DegreesOfFreedom<World> const world_degrees_of_freedom =
        from_navigation_frame_to_world_at_current_time(
            navigation_degrees_of_freedom);
    trajectory->Append(it.time(), world_degrees_of_freedom);
  }
  VLOG(1) << "Returning a " << trajectory->Size() << "-point trajectory";
  return trajectory;
}

AffineMap<Barycentric, World, Length, OrthogonalMap>
Renderer::BarycentricToWorld(
    Instant const& time,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return AffineMap<Barycentric, World, Length, OrthogonalMap>(
      sun_->current_position(time),
      sun_world_position,
      BarycentricToWorld(planetarium_rotation));
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

AffineMap<World, Barycentric, Length, OrthogonalMap>
Renderer::WorldToBarycentric(
    Instant const& time,
    Position<World> const& sun_world_position,
    Rotation<Barycentric, AliceSun> const& planetarium_rotation) const {
  return AffineMap<World, Barycentric, Length, OrthogonalMap>(
      sun_world_position,
      sun_->current_position(time),
      WorldToBarycentric(planetarium_rotation));
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
