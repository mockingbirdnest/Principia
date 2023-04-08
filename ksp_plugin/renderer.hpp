#pragma once

#include <functional>
#include <memory>
#include <optional>

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _renderer {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_affine_map;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_quantities;

class Renderer {
 public:
  Renderer(not_null<Celestial const*> sun,
           not_null<std::unique_ptr<NavigationFrame>> plotting_frame);

  virtual ~Renderer() = default;

  // Changes the plotting frame of the renderer.
  virtual void SetPlottingFrame(
      not_null<std::unique_ptr<NavigationFrame>> plotting_frame);

  // Returns the current plotting frame.  This may not be the last set by
  // |SetPlottingFrame| if it is overridden by a target vessel.
  virtual not_null<NavigationFrame const*> GetPlottingFrame() const;

  // Overrides the current plotting frame with one that is centred on the given
  // |vessel|.
  virtual void SetTargetVessel(
      not_null<Vessel*> vessel,
      not_null<Celestial const*> celestial,
      not_null<Ephemeris<Barycentric> const*> ephemeris);

  // Reverts to frame last set by |SetPlottingFrame|.  The second version only
  // has an effect if the given |vessel| is the current target vessel.
  virtual void ClearTargetVessel();
  virtual void ClearTargetVesselIf(not_null<Vessel*> vessel);

  // Determines if there is a target vessel and returns it.
  virtual bool HasTargetVessel() const;
  virtual Vessel& GetTargetVessel();
  virtual Vessel const& GetTargetVessel() const;

  // Returns a trajectory in |World| corresponding to the trajectory defined by
  // |begin| and |end|, as seen in the current plotting frame.  In this function
  // and others in this class, |sun_world_position| is the current position of
  // the sun in |World| space as returned by |Planetarium.fetch.Sun.position|;
  // it is used to define the relation between |WorldSun| and |World|.
  virtual DiscreteTrajectory<World>
  RenderBarycentricTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Barycentric>::iterator const& begin,
      DiscreteTrajectory<Barycentric>::iterator const& end,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Returns a trajectory in the current plotting frame corresponding to the
  // trajectory defined by |begin| and |end|.  If there is a target vessel, its
  // prediction must not be empty.
  virtual DiscreteTrajectory<Navigation>
  RenderBarycentricTrajectoryInPlotting(
      DiscreteTrajectory<Barycentric>::iterator const& begin,
      DiscreteTrajectory<Barycentric>::iterator const& end) const;

  // Returns a trajectory in |World| corresponding to the trajectory defined by
  // |begin| and |end| in the current plotting frame.
  virtual DiscreteTrajectory<World>
  RenderPlottingTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Navigation>::iterator const& begin,
      DiscreteTrajectory<Navigation>::iterator const& end,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Coordinate transforms.

  virtual RigidMotion<Barycentric, Navigation> BarycentricToPlotting(
      Instant const& time) const;

  virtual RigidTransformation<Barycentric, World> BarycentricToWorld(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Barycentric, World> BarycentricToWorld(
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Barycentric, WorldSun> BarycentricToWorldSun(
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Converts from the Frenet frame of the manœuvre's initial time in the
  // plotted frame to the |World| coordinates.
  virtual OrthogonalMap<Frenet<Navigation>, World> FrenetToWorld(
      Instant const& time,
      NavigationManœuvre const& manœuvre,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Converts from the Frenet frame of the vessel's free-falling trajectory in
  // the plotted frame to the |World| coordinates.
  virtual OrthogonalMap<Frenet<Navigation>, World> FrenetToWorld(
      Vessel const& vessel,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Converts from the Frenet frame of the vessel's free-falling trajectory in
  // the given |navigation_frame| to the |World| coordinates.
  virtual OrthogonalMap<Frenet<Navigation>, World> FrenetToWorld(
      Vessel const& vessel,
      NavigationFrame const& navigation_frame,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Navigation, Barycentric> PlottingToBarycentric(
      Instant const& time) const;

  virtual RigidTransformation<Navigation, World> PlottingToWorld(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Navigation, World> PlottingToWorld(
      Instant const& time,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual RigidTransformation<World, Barycentric> WorldToBarycentric(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<World, Barycentric> WorldToBarycentric(
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual RigidTransformation<World, Navigation> WorldToPlotting(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual Rotation<CameraCompensatedReference, World> CameraReferenceRotation(
      Instant const& time,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation,
      Rotation<CameraCompensatedReference, CameraReference> const&
          camera_compensation) const;

  virtual void WriteToMessage(not_null<serialization::Renderer*> message) const;

  static not_null<std::unique_ptr<Renderer>> ReadFromMessage(
      serialization::Renderer const& message,
      not_null<Celestial const*> sun,
      not_null<Ephemeris<Barycentric> const*> ephemeris);

 private:
  struct Target {
    Target(not_null<Vessel*> vessel,
           not_null<Celestial const*> celestial,
           not_null<Ephemeris<Barycentric> const*> ephemeris);
    not_null<Vessel*> const vessel;
    not_null<Celestial const*> const celestial;
    not_null<std::unique_ptr<NavigationFrame>> const target_frame;
  };

  not_null<Celestial const*> const sun_;

  not_null<std::unique_ptr<NavigationFrame>> plotting_frame_;

  std::optional<Target> target_;
};

}  // namespace internal

using internal::Renderer;

}  // namespace _renderer
}  // namespace ksp_plugin
}  // namespace principia

namespace principia::ksp_plugin {
using namespace principia::ksp_plugin::_renderer;
}  // namespace principia::ksp_plugin
