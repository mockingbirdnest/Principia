#include "ksp_plugin/interface.hpp"

#include <algorithm>
#include <limits>

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/rp2_point.hpp"
#include "geometry/space_transformations.hpp"
#include "glog/logging.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/renderer.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::geometry::_affine_map;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_perspective;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_rp2_point;
using namespace principia::geometry::_space_transformations;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_iterators;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

Planetarium* __cdecl principia__PlanetariumCreate(
    Plugin const* const plugin,
    XYZ const sun_world_position,
    XYZ const xyz_opengl_camera_x_in_world,
    XYZ const xyz_opengl_camera_y_in_world,
    XYZ const xyz_opengl_camera_z_in_world,
    XYZ const xyz_camera_position_in_world,
    double const focal,
    double const field_of_view,
    double const inverse_scale_factor,
    XYZ const scaled_space_origin) {
  journal::Method<journal::PlanetariumCreate> m({plugin,
                                                 sun_world_position,
                                                 xyz_opengl_camera_x_in_world,
                                                 xyz_opengl_camera_y_in_world,
                                                 xyz_opengl_camera_z_in_world,
                                                 xyz_camera_position_in_world,
                                                 focal,
                                                 field_of_view,
                                                 inverse_scale_factor,
                                                 scaled_space_origin});
  Renderer const& renderer = CHECK_NOTNULL(plugin)->renderer();

  Multivector<double, World, 1> const opengl_camera_x_in_world(
      FromXYZ(xyz_opengl_camera_x_in_world));
  Multivector<double, World, 1> const opengl_camera_y_in_world(
      FromXYZ(xyz_opengl_camera_y_in_world));
  Multivector<double, World, 2> const opengl_camera_z_in_world(
      FromXYZ(xyz_opengl_camera_z_in_world));
  // Note the minus sign for z below because our convention with respect to the
  // orientation of z is opposite that of OpenGL.
  Rotation<Camera, World> const camera_to_world_rotation(
      opengl_camera_x_in_world,
      opengl_camera_y_in_world,
      -opengl_camera_z_in_world);
  Position<World> const camera_position_in_world =
      FromXYZ<Position<World>>(xyz_camera_position_in_world);

  RigidTransformation<Camera, World> const camera_to_world_affine_map(
      Camera::origin,
      camera_position_in_world,
      camera_to_world_rotation.Forget<OrthogonalMap>());
  Similarity<World, Navigation> const
      world_to_plotting_affine_map =
          renderer.WorldToPlotting(plugin->CurrentTime(),
                                   FromXYZ<Position<World>>(sun_world_position),
                                   plugin->PlanetariumRotation());

  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1.0,
      /*angular_resolution=*/0.4 * ArcMinute,
      field_of_view * Radian);
  Perspective<Navigation, Camera> perspective(
      world_to_plotting_affine_map *
          camera_to_world_affine_map.Forget<Similarity>(),
      focal * Metre);

  auto const plotting_to_scaled_space =
      [plotting_to_world = world_to_plotting_affine_map.Inverse(),
       scaled_space_origin = FromXYZ<Position<World>>(scaled_space_origin),
       inverse_scale_factor = inverse_scale_factor * (1 / Metre)](
          Position<Navigation> const& plotted_point) {
        return ScaledSpacePoint::FromCoordinates(
            ((plotting_to_world(plotted_point) - scaled_space_origin) *
             inverse_scale_factor).coordinates());
      };
  return m.Return(
      plugin->NewPlanetarium(
          parameters,
          perspective,
          plotting_to_scaled_space).release());
}

void __cdecl principia__PlanetariumDelete(
    Planetarium const** const planetarium) {
  journal::Method<journal::PlanetariumDelete> m({planetarium}, {planetarium});
  CHECK_NOTNULL(planetarium);
  TakeOwnership(planetarium);
  return m.Return();
}

// Fills the array of size |vertices_size| at |vertices| with vertices for the
// rendering of the segment with the given index in the flight plan of the
// vessel with the given GUID.
void __cdecl principia__PlanetariumPlotFlightPlanSegment(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index,
    ScaledSpacePoint* const vertices,
    int const vertices_size,
    int* const vertex_count) {
  journal::Method<journal::PlanetariumPlotFlightPlanSegment> m(
      {planetarium, plugin, vessel_guid, index, vertices, vertices_size},
      {vertex_count});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  *vertex_count = 0;

  Vessel const& vessel = *plugin->GetVessel(vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  auto const segment = vessel.flight_plan().GetSegment(index);
  // TODO(egg): this is ugly; we should centralize rendering.
  // If this is a burn and we cannot render the beginning of the burn, we
  // render none of it, otherwise we try to render the Frenet trihedron at the
  // start and we fail.
  if (index % 2 == 0 ||
      segment->empty() ||
      segment->front().time >= plugin->renderer().GetPlottingFrame()->t_min()) {
    planetarium->PlotMethod3(
        *segment, segment->begin(), segment->end(),
        plugin->CurrentTime(),
        /*reverse=*/false,
        [vertices, vertex_count](ScaledSpacePoint const& vertex) {
          vertices[(*vertex_count)++] = vertex;
        },
        vertices_size);
  }
  return m.Return();
}

// Fills the array of size |vertices_size| at |vertices| with vertices for the
// rendered prediction of the vessel with the given GUID.
void __cdecl principia__PlanetariumPlotPrediction(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    char const* const vessel_guid,
    ScaledSpacePoint* const vertices,
    int const vertices_size,
    int* const vertex_count) {
  journal::Method<journal::PlanetariumPlotPrediction> m(
      {planetarium, plugin, vessel_guid, vertices, vertices_size},
      {vertex_count});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  *vertex_count = 0;

  auto const prediction = plugin->GetVessel(vessel_guid)->prediction();
  planetarium->PlotMethod3(
      *prediction, prediction->begin(), prediction->end(),
      plugin->CurrentTime(),
      /*reverse=*/false,
      [vertices, vertex_count](ScaledSpacePoint const& vertex) {
        vertices[(*vertex_count)++] = vertex;
      },
      vertices_size);
  return m.Return();
}

// Fills the array of size |vertices_size| at |vertices| with vertices for the
// rendered past trajectory of the vessel with the given GUID; the
// trajectory goes back |max_history_length| seconds before the present time (or
// to the earliest time available if the relevant |t_min| is more recent).
void __cdecl principia__PlanetariumPlotPsychohistory(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const max_history_length,
    ScaledSpacePoint* const vertices,
    int const vertices_size,
    int* const vertex_count) {
  journal::Method<journal::PlanetariumPlotPsychohistory> m(
      {planetarium,
       plugin,
       vessel_guid,
       max_history_length,
       vertices,
       vertices_size},
      {vertex_count});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  *vertex_count = 0;

  // Do not plot the psychohistory when there is a target vessel as it is
  // misleading.
  if (plugin->renderer().HasTargetVessel()) {
    return m.Return();
  } else {
    auto const vessel = plugin->GetVessel(vessel_guid);
    auto const& trajectory = vessel->trajectory();
    auto const& psychohistory = vessel->psychohistory();

    Instant const desired_first_time =
        plugin->CurrentTime() - max_history_length * Second;

    // Since we would want to plot starting from |desired_first_time|, ask the
    // reanimator to reconstruct the past.  That may take a while, during which
    // time the history will be shorter than desired.
    vessel->RequestReanimation(desired_first_time);

    planetarium->PlotMethod3(
        trajectory,
        trajectory.lower_bound(desired_first_time),
        psychohistory->end(),
        /*now=*/plugin->CurrentTime(),
        /*reverse=*/true,
        [vertices, vertex_count](ScaledSpacePoint const& vertex) {
          vertices[(*vertex_count)++] = vertex;
        },
        vertices_size);
    return m.Return();
  }
}

// Fills the array of size |vertices_size| at |vertices| with vertices for the
// rendered past trajectory of the celestial with the given index; the
// trajectory goes back |max_history_length| seconds before the present time (or
// to the earliest time available if the relevant |t_min| is more recent).
void __cdecl principia__PlanetariumPlotCelestialPastTrajectory(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const celestial_index,
    double const max_history_length,
    ScaledSpacePoint* const vertices,
    int const vertices_size,
    double* const minimal_distance_from_camera,
    int* const vertex_count) {
  journal::Method<journal::PlanetariumPlotCelestialPastTrajectory> m(
      {planetarium,
       plugin,
       celestial_index,
       max_history_length,
       vertices,
       vertices_size},
      {minimal_distance_from_camera, vertex_count});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  *vertex_count = 0;

  // Do not plot the past when there is a target vessel as it is misleading.
  if (plugin->renderer().HasTargetVessel()) {
    *minimal_distance_from_camera = std::numeric_limits<double>::infinity();
    return m.Return();
  } else {
    auto const& celestial_trajectory =
        plugin->GetCelestial(celestial_index).trajectory();
    Instant const desired_first_time =
        plugin->CurrentTime() - max_history_length * Second;

    // Since we would want to plot starting from |desired_first_time|, ask the
    // reanimator to reconstruct the past.  That may take a while, during which
    // time the history will be shorter than desired.
    plugin->RequestReanimation(desired_first_time);

    Instant const first_time =
        std::max(desired_first_time, celestial_trajectory.t_min());
    Length minimal_distance;
    planetarium->PlotMethod3(
        celestial_trajectory,
        first_time,
        /*last_time=*/plugin->CurrentTime(),
        /*now=*/plugin->CurrentTime(),
        /*reverse=*/true,
        [vertices, vertex_count](ScaledSpacePoint const& vertex) {
          vertices[(*vertex_count)++] = vertex;
        },
        vertices_size,
        &minimal_distance);
    *minimal_distance_from_camera = minimal_distance / Metre;
    return m.Return();
  }
}

// Fills the array of size |vertices_size| at |vertices| with vertices for the
// rendered future trajectory of the celestial with the given index; the
// trajectory goes as far as the furthest of the final time of the prediction or
// that of the flight plan.
void __cdecl principia__PlanetariumPlotCelestialFutureTrajectory(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const celestial_index,
    char const* const vessel_guid,
    ScaledSpacePoint* const vertices,
    int const vertices_size,
    double* const minimal_distance_from_camera,
    int* const vertex_count) {
  journal::Method<journal::PlanetariumPlotCelestialFutureTrajectory> m(
      {planetarium,
       plugin,
       celestial_index,
       vessel_guid,
       vertices,
       vertices_size},
      {minimal_distance_from_camera, vertex_count});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  *vertex_count = 0;

  // Do not plot the past when there is a target vessel as it is misleading.
  // TODO(egg): This is the future, not the past!
  if (plugin->renderer().HasTargetVessel()) {
    *minimal_distance_from_camera = std::numeric_limits<double>::infinity();
    return m.Return();
  } else {
    auto const& vessel = *plugin->GetVessel(vessel_guid);
    Instant const prediction_final_time = vessel.prediction()->t_max();
    Instant const final_time =
        vessel.has_flight_plan()
            ? std::max(vessel.flight_plan().actual_final_time(),
                       prediction_final_time)
            : prediction_final_time;
    auto const& celestial_trajectory =
        plugin->GetCelestial(celestial_index).trajectory();
    // No need to request reanimation here because the current time of the
    // plugin is necessarily covered.
    Length minimal_distance;
    planetarium->PlotMethod3(
        celestial_trajectory,
        /*first_time=*/plugin->CurrentTime(),
        /*last_time=*/final_time,
        /*now=*/plugin->CurrentTime(),
        /*reverse=*/false,
        [vertices, vertex_count](ScaledSpacePoint const& vertex) {
          vertices[(*vertex_count)++] = vertex;
        },
        vertices_size,
        &minimal_distance);
    *minimal_distance_from_camera = minimal_distance / Metre;
    return m.Return();
  }
}

}  // namespace interface
}  // namespace principia
