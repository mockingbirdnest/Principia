#pragma once

#include <experimental/optional>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using base::not_null;
using geometry::AffineMap;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Rotation;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::Frenet;
using physics::RigidMotion;
using quantities::Length;

class Renderer {
 public:
  Renderer(not_null<Celestial const*> sun,
           not_null<std::unique_ptr<NavigationFrame>> plotting_frame);

  // Changes the plotting frame of the renderer.
  virtual void SetPlottingFrame(
      not_null<std::unique_ptr<NavigationFrame>> plotting_frame);

  // Returns the current plotting frame.  This may not be the last set by
  // |SetPlottingFrame| if it was overridden by a target vessel.
  virtual not_null<NavigationFrame const*> GetPlottingFrame() const;

  // Overrides the current plotting frame with one that is centred on the given
  // |vessel|.  When using the operations below with a target vessel, the client
  // must ensure that the |time| is covered by the vessel's prediction.
  virtual void SetTargetVessel(
      not_null<Vessel*> vessel,
      not_null<Celestial const*> celestial,
      not_null<Ephemeris<Barycentric> const*> const ephemeris);

  // Reverts to the previously active plotting frame.
  virtual void ClearTargetVessel();

  // Determines if there is a target vessel and returns it.
  virtual bool HasTargetVessel() const;
  virtual Vessel const& GetTargetVessel() const;

  // Returns a trajectory in |World| corresponding to the trajectory defined by
  // |begin| and |end|, as seen in the current plotting frame.  In this function
  // and others in this class, |sun_world_position| is the current position of
  // the sun in |World| space as returned by |Planetarium.fetch.Sun.position|;
  // it is used to define the relation between |WorldSun| and |World|.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderBarycentricTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Returns a trajectory in the current plotting frame corresponding to the
  // trajectory defined by |begin| and |end|.  If there is a target vessel, its
  // prediction must cover the time corresponding to |end|.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<Navigation>>>
  RenderBarycentricTrajectoryInPlotting(
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end) const;

  // Returns a trajectory in |World| corresponding to the trajectory defined by
  // |begin| and |end| in the current plotting frame.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderPlottingTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Navigation>::Iterator const& begin,
      DiscreteTrajectory<Navigation>::Iterator const& end,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Coordinate transforms.

  virtual RigidMotion<Barycentric, Navigation> BarycentricToPlotting(
      Instant const& time) const;

  virtual RigidMotion<Barycentric, World> BarycentricToWorld(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Barycentric, World> BarycentricToWorld(
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Barycentric, WorldSun> BarycentricToWorldSun(
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  // Converts from the Frenet frame of the vessel's free-falling trajectory in
  // the plotted frame to the |World| coordinates.
  virtual OrthogonalMap<Frenet<Navigation>, World> FrenetToWorld(
      Vessel const& vessel,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Navigation, Barycentric> PlottingToBarycentric(
      Instant const& time) const;

  virtual RigidMotion<Navigation, World> PlottingToWorld(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<Navigation, World> PlottingToWorld(
      Instant const& time,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual RigidMotion<World, Barycentric> WorldToBarycentric(
      Instant const& time,
      Position<World> const& sun_world_position,
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

  virtual OrthogonalMap<World, Barycentric> WorldToBarycentric(
      Rotation<Barycentric, AliceSun> const& planetarium_rotation) const;

 private:
  not_null<Celestial const*> const sun_;

  not_null<std::unique_ptr<NavigationFrame>> plotting_frame_;

  struct Target {
    Target(not_null<Vessel*> vessel,
           not_null<Celestial const*> celestial,
           not_null<Ephemeris<Barycentric> const*> ephemeris);
    not_null<Vessel*> const vessel;
    not_null<Celestial const*> const celestial;
    not_null<std::unique_ptr<NavigationFrame>> const target_frame;
  };
  std::experimental::optional<Target> target_;
};

}  // namespace internal_renderer
}  // namespace ksp_plugin
}  // namespace principia
