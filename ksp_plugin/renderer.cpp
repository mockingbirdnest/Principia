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
  if (!target_ || target_->vessel != vessel ||
      target_->celestial != celestial) {
    target_.emplace(vessel, celestial, ephemeris);
  }
  // Make sure that the current time is covered by the prediction.
  if (current_time_ > target_->vessel->prediction().t_max()) {
    target_->vessel->UpdatePrediction(current_time_ + prediction_length_);
  }
}

void Renderer::ClearTargetVessel() {
  target_ = std::experimental::nullopt;
}

not_null<std::unique_ptr<DiscreteTrajectory<World>>>
Renderer::RenderBarycentricTrajectoryInWorld(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  auto const trajectory_in_navigation =
      RenderBarycentricTrajectoryInNavigation(begin, end);
  auto trajectory_in_world =
      RenderNavigationTrajectoryInWorld(trajectory_in_navigation->Begin(),
                                        trajectory_in_navigation->End(),
                                        sun_world_position);
  return trajectory_in_world;
}

void Renderer::ComputeAndRenderApsides(
    Celestial const& celestial,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    std::unique_ptr<DiscreteTrajectory<World>>& apoapsides,
    std::unique_ptr<DiscreteTrajectory<World>>& periapsides) const {
  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(celestial.trajectory(),
                 begin,
                 end,
                 apoapsides_trajectory,
                 periapsides_trajectory);
  apoapsides =
      RenderBarycentricTrajectoryInWorld(apoapsides_trajectory.Begin(),
                                         apoapsides_trajectory.End(),
                                         sun_world_position);
  periapsides =
      RenderBarycentricTrajectoryInWorld(periapsides_trajectory.Begin(),
                                         periapsides_trajectory.End(),
                                         sun_world_position);
}

void Renderer::ComputeAndRenderClosestApproaches(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    std::unique_ptr<DiscreteTrajectory<World>>& closest_approaches) const {
  CHECK(target_);

  DiscreteTrajectory<Barycentric> apoapsides_trajectory;
  DiscreteTrajectory<Barycentric> periapsides_trajectory;
  ComputeApsides(target_->vessel->prediction(),
                 begin,
                 end,
                 apoapsides_trajectory,
                 periapsides_trajectory);
  closest_approaches =
      RenderBarycentricTrajectoryInWorld(periapsides_trajectory.Begin(),
                                         periapsides_trajectory.End(),
                                         sun_world_position);
}

void Renderer::ComputeAndRenderNodes(
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Position<World> const& sun_world_position,
    std::unique_ptr<DiscreteTrajectory<World>>& ascending,
    std::unique_ptr<DiscreteTrajectory<World>>& descending) const {
  CHECK(target_);
  auto const trajectory_in_navigation =
      RenderBarycentricTrajectoryInNavigation(begin, end);
  DiscreteTrajectory<Navigation> ascending_trajectory;
  DiscreteTrajectory<Navigation> descending_trajectory;
  // The so-called North is orthogonal to the plane of the trajectory.
  ComputeNodes(trajectory_in_navigation->Begin(),
               trajectory_in_navigation->End(),
               Vector<double, Navigation>({0, 0, 1}),
               ascending_trajectory,
               descending_trajectory);
  ascending = RenderNavigationTrajectoryInWorld(ascending_trajectory.Begin(),
                                                ascending_trajectory.End(),
                                                sun_world_position);
  descending = RenderNavigationTrajectoryInWorld(descending_trajectory.Begin(),
                                                 descending_trajectory.End(),
                                                 sun_world_position);
}

AffineMap<Barycentric, World, Length, OrthogonalMap>
Renderer::BarycentricToWorld(Position<World> const& sun_world_position) const {
  return AffineMap<Barycentric, World, Length, OrthogonalMap>(
      sun_->current_position(current_time_),
      sun_world_position,
      BarycentricToWorld());
}

OrthogonalMap<Barycentric, World> Renderer::BarycentricToWorld() const {
  return OrthogonalMap<WorldSun, World>::Identity() * BarycentricToWorldSun();
}

OrthogonalMap<Barycentric, WorldSun> Renderer::BarycentricToWorldSun() const {
  return sun_looking_glass.Inverse().Forget() * PlanetariumRotation().Forget();
}

AffineMap<World, Barycentric, Length, OrthogonalMap>
Renderer::WorldToBarycentric(Position<World> const& sun_world_position) const {
  return AffineMap<World, Barycentric, Length, OrthogonalMap>(
      sun_world_position,
      sun_->current_position(current_time_),
      WorldToBarycentric());
}

OrthogonalMap<World, Barycentric> Renderer::WorldToBarycentric() const {
  return BarycentricToWorld().Inverse();
}

not_null<std::unique_ptr<DiscreteTrajectory<Navigation>>>
Renderer::RenderBarycentricTrajectoryInNavigation(
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
    target_->vessel->UpdatePrediction(current_time_ + prediction_length_);
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
    DiscreteTrajectory<Navigation>::Iterator const& begin,
    DiscreteTrajectory<Navigation>::Iterator const& end,
    Position<World> const& sun_world_position) const {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<World>>();

  NavigationFrame const& plotting_frame = *GetPlottingFrame();

  RigidMotion<Navigation, World> from_navigation_frame_to_world_at_current_time(
      /*rigid_transformation=*/BarycentricToWorld(sun_world_position) *
          plotting_frame.FromThisFrameAtTime(current_time_)
              .rigid_transformation(),
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
