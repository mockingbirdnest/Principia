#pragma once

#include <experimental/optional>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_renderer {

using base::not_null;
using geometry::AffineMap;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using quantities::Length;

class Renderer {
 public:
  Renderer(not_null<Celestial const*> sun,
           not_null<std::unique_ptr<NavigationFrame>> plotting_frame);

  virtual void SetPlottingFrame(
      not_null<std::unique_ptr<NavigationFrame>> plotting_frame);
  virtual not_null<NavigationFrame const*> GetPlottingFrame() const;

  virtual void SetTargetVessel(
      Instant const& time,
      not_null<Vessel*> vessel,
      not_null<Celestial const*> celestial,
      not_null<Ephemeris<Barycentric> const*> const ephemeris);
  virtual void ClearTargetVessel();
  virtual bool HasTargetVessel() const;
  virtual Vessel const& GetTargetVessel() const;

  // Returns a |Trajectory| object corresponding to the trajectory defined by
  // |begin| and |end|, as seen in the current |plotting_frame_|.
  virtual not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderBarycentricTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Position<World> const& sun_world_position) const;

  // Converts a trajectory from |Barycentric| to |Navigation|.
  not_null<std::unique_ptr<DiscreteTrajectory<Navigation>>>
  RenderBarycentricTrajectoryInNavigation(
      Instant const& time,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end) const;

  // Converts a trajectory from |Navigation| to |World|.  |sun_world_position|
  // is the current position of the sun in |World| space as returned by
  // |Planetarium.fetch.Sun.position|.  It is used to define the relation
  // between |WorldSun| and |World|.
  not_null<std::unique_ptr<DiscreteTrajectory<World>>>
  RenderNavigationTrajectoryInWorld(
      Instant const& time,
      DiscreteTrajectory<Navigation>::Iterator const& begin,
      DiscreteTrajectory<Navigation>::Iterator const& end,
      Position<World> const& sun_world_position) const;

  // Coordinate transforms.
  virtual AffineMap<Barycentric, World, Length, OrthogonalMap>
  BarycentricToWorld(Instant const& time,
                     Position<World> const& sun_world_position) const;
  virtual OrthogonalMap<Barycentric, World> BarycentricToWorld() const;
  virtual OrthogonalMap<Barycentric, WorldSun> BarycentricToWorldSun() const;
  virtual AffineMap<World, Barycentric, Length, OrthogonalMap>
  WorldToBarycentric(Instant const& time,
                     Position<World> const& sun_world_position) const;
  virtual OrthogonalMap<World, Barycentric> WorldToBarycentric() const;

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