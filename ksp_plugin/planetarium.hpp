
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rp2_point.hpp"
#include "geometry/sphere.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Perspective;
using geometry::RP2Line;
using geometry::RP2Point;
using geometry::Segment;
using geometry::Sphere;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::RigidMotion;
using quantities::Length;

// A planetarium is an ephemeris together with a perspective.  In this setting
// it is possible to draw trajectories in the projective plane.
class Planetarium final {
 public:
  // TODO(phl): All this Navigation is weird.  Should it be named Plotting?
  // In particular Navigation vs. NavigationFrame is a mess.
  Planetarium(Perspective<Navigation, Camera, Length, OrthogonalMap> const&
                  perspective,
              not_null<Ephemeris<Barycentric> const*> ephemeris,
              not_null<NavigationFrame*> plotting_frame);

  // A no-op method that just returns all the points in the trajectory defined
  // by |begin| and |end|.
  std::vector<RP2Line<Length, Camera>> PlotMethod0(
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Instant const& now) const;

  // A naïve method that doesn't pay any attention to the perspective but tries
  // to ensure that the points before the perspective are separated by less than
  // |tolerance|.
  std::vector<RP2Line<Length, Camera>> PlotMethod1(
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end,
      Instant const& now,
      Length const& tolerance) const;

 private:
  // Computes the coordinates of the spheres that represent the |ephemeris_|
  // bodies.  These coordinates are in the |plotting_frame_| at time |now|.
  std::vector<Sphere<Length, Navigation>> ComputePlottableSpheres(
      Instant const& now) const;

  //TODO(phl):comment
  std::vector<Segment<Displacement<Navigation>>> ComputePlottableSegments(
      const std::vector<Sphere<Length, Navigation>>& plottable_spheres,
      DiscreteTrajectory<Barycentric>::Iterator const& begin,
      DiscreteTrajectory<Barycentric>::Iterator const& end) const;

  Perspective<Navigation, Camera, Length, OrthogonalMap> const
      perspective_;
  not_null<Ephemeris<Barycentric> const*> const ephemeris_;
  not_null<NavigationFrame*> const plotting_frame_;
};

}  // namespace internal_planetarium

using internal_planetarium::Planetarium;

}  // namespace ksp_plugin
}  // namespace principia
