
#include "ksp_plugin/interface.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/rp2_point.hpp"
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

using geometry::AffineMap;
using geometry::Multivector;
using geometry::OrthogonalMap;
using geometry::Perspective;
using geometry::RigidTransformation;
using geometry::Rotation;
using geometry::RP2Lines;
using ksp_plugin::Camera;
using ksp_plugin::Navigation;
using ksp_plugin::Planetarium;
using ksp_plugin::Renderer;
using ksp_plugin::TypedIterator;
using physics::DiscreteTrajectory;
using quantities::Length;
using quantities::si::ArcMinute;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Radian;

namespace {

RP2Lines<Length, Camera> PlotMethodN(
    Planetarium const& planetarium,
    int const method,
    DiscreteTrajectory<Barycentric>::Iterator const& begin,
    DiscreteTrajectory<Barycentric>::Iterator const& end,
    Instant const& now,
    bool const reverse) {
  switch (method) {
    case 0:
      return planetarium.PlotMethod0(begin, end, now, reverse);
    case 1:
      return planetarium.PlotMethod1(begin, end, now, reverse);
    case 2:
      return planetarium.PlotMethod2(begin, end, now, reverse);
    default:
      LOG(FATAL) << "Unexpected method " << method;
      base::noreturn();
  }
}

}  // namespace

Planetarium* principia__PlanetariumCreate(
    Plugin const* const plugin,
    XYZ const sun_world_position,
    XYZ const xyz_opengl_camera_x_in_world,
    XYZ const xyz_opengl_camera_y_in_world,
    XYZ const xyz_opengl_camera_z_in_world,
    XYZ const xyz_camera_position_in_world,
    double const focal,
    double const field_of_view) {
  journal::Method<journal::PlanetariumCreate> m({plugin,
                                                 sun_world_position,
                                                 xyz_opengl_camera_x_in_world,
                                                 xyz_opengl_camera_y_in_world,
                                                 xyz_opengl_camera_z_in_world,
                                                 xyz_camera_position_in_world,
                                                 focal,
                                                 field_of_view});
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

  RigidTransformation<Camera, World> const
      camera_to_world_affine_map(Camera::origin,
                                 camera_position_in_world,
                                 camera_to_world_rotation.Forget());
  RigidTransformation<World, Navigation> const
      world_to_plotting_affine_map =
          renderer.WorldToPlotting(plugin->CurrentTime(),
                                   FromXYZ<Position<World>>(sun_world_position),
                                   plugin->PlanetariumRotation());

  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1.0,
      /*angular_resolution=*/0.4 * ArcMinute,
      field_of_view * Radian);
  Perspective<Navigation, Camera> perspective(
      world_to_plotting_affine_map * camera_to_world_affine_map,
      focal * Metre);

  return m.Return(plugin->NewPlanetarium(parameters, perspective).release());
}

void principia__PlanetariumDelete(
    Planetarium const** const planetarium) {
  journal::Method<journal::PlanetariumDelete> m({planetarium}, {planetarium});
  CHECK_NOTNULL(planetarium);
  TakeOwnership(planetarium);
  return m.Return();
}

Iterator* principia__PlanetariumPlotFlightPlanSegment(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const method,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::PlanetariumPlotFlightPlanSegment> m({planetarium,
                                                                plugin,
                                                                method,
                                                                vessel_guid,
                                                                index});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  Vessel const& vessel = *plugin->GetVessel(vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  DiscreteTrajectory<Barycentric>::Iterator segment_begin;
  DiscreteTrajectory<Barycentric>::Iterator segment_end;
  vessel.flight_plan().GetSegment(index, segment_begin, segment_end);
  RP2Lines<Length, Camera> rp2_lines;
  // TODO(egg): this is ugly; we should centralize rendering.
  // If this is a burn and we cannot render the beginning of the burn, we
  // render none of it, otherwise we try to render the Frenet trihedron at the
  // start and we fail.
  if (index % 2 == 0 ||
      segment_begin == segment_end ||
      segment_begin->time >= plugin->renderer().GetPlottingFrame()->t_min()) {
     rp2_lines = PlotMethodN(*planetarium,
                             method,
                             segment_begin,
                             segment_end,
                             plugin->CurrentTime(),
                             /*reverse=*/false);
  }
  return m.Return(new TypedIterator<RP2Lines<Length, Camera>>(rp2_lines));
}

Iterator* principia__PlanetariumPlotPrediction(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const method,
    char const* const vessel_guid) {
  journal::Method<journal::PlanetariumPlotPrediction> m({planetarium,
                                                         plugin,
                                                         method,
                                                         vessel_guid});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);
  auto const& prediction = plugin->GetVessel(vessel_guid)->prediction();
  auto const rp2_lines = PlotMethodN(*planetarium,
                                     method,
                                     prediction.Fork(),
                                     prediction.end(),
                                     plugin->CurrentTime(),
                                     /*reverse=*/false);
  return m.Return(new TypedIterator<RP2Lines<Length, Camera>>(rp2_lines));
}

Iterator* principia__PlanetariumPlotPsychohistory(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const method,
    char const* const vessel_guid,
    double const max_history_length) {
  journal::Method<journal::PlanetariumPlotPsychohistory> m(
      {planetarium, plugin, method, vessel_guid, max_history_length});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);

  // Do not plot the psychohistory when there is a target vessel as it is
  // misleading.
  if (plugin->renderer().HasTargetVessel()) {
    return m.Return(new TypedIterator<RP2Lines<Length, Camera>>({}));
  } else {
    auto const& psychohistory = plugin->GetVessel(vessel_guid)->psychohistory();
    auto const rp2_lines =
        PlotMethodN(*planetarium,
                    method,
                    psychohistory.LowerBound(plugin->CurrentTime() -
                                             max_history_length * Second),
                    psychohistory.end(),
                    plugin->CurrentTime(),
                    /*reverse=*/true);
    return m.Return(new TypedIterator<RP2Lines<Length, Camera>>(rp2_lines));
  }
}

Iterator* principia__PlanetariumPlotCelestialTrajectoryForPsychohistory(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const celestial_index,
    char const* const vessel_guid,
    double const max_history_length) {
  journal::Method<journal::PlanetariumPlotCelestialTrajectoryForPsychohistory>
      m({planetarium, plugin, celestial_index, vessel_guid, max_history_length});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);

  // Do not plot the past when there is a target vessel as it is misleading.
  if (plugin->renderer().HasTargetVessel()) {
    return m.Return(new TypedIterator<RP2Lines<Length, Camera>>({}));
  } else {
    auto const& celestial_trajectory =
        plugin->GetCelestial(celestial_index).trajectory();
    Instant const first_time = std::max(
        plugin->CurrentTime() - max_history_length * Second,
        vessel_guid == nullptr
            ? celestial_trajectory.t_min()
            : plugin->GetVessel(vessel_guid)->psychohistory().front().time);
    auto const rp2_lines =
        planetarium->PlotMethod2(celestial_trajectory,
                                 first_time,
                                 /*last_time=*/plugin->CurrentTime(),
                                 /*now=*/plugin->CurrentTime(),
                                 /*reverse=*/true);
    return m.Return(new TypedIterator<RP2Lines<Length, Camera>>(rp2_lines));
  }
}

Iterator* principia__PlanetariumPlotCelestialTrajectoryForPredictionOrFlightPlan(
    Planetarium const* const planetarium,
    Plugin const* const plugin,
    int const celestial_index,
    char const* const vessel_guid) {
  journal::Method<
      journal::PlanetariumPlotCelestialTrajectoryForPredictionOrFlightPlan>
      m({planetarium, plugin, celestial_index, vessel_guid});
  CHECK_NOTNULL(plugin);
  CHECK_NOTNULL(planetarium);

  // Do not plot the past when there is a target vessel as it is misleading.
  if (plugin->renderer().HasTargetVessel()) {
    return m.Return(new TypedIterator<RP2Lines<Length, Camera>>({}));
  } else {
    auto const& vessel = *plugin->GetVessel(vessel_guid);
    Instant const prediction_final_time = vessel.prediction().back().time;
    Instant const final_time =
        vessel.has_flight_plan()
            ? std::max(vessel.flight_plan().actual_final_time(),
                       prediction_final_time)
            : prediction_final_time;
    auto const& celestial_trajectory =
        plugin->GetCelestial(celestial_index).trajectory();
    auto const rp2_lines =
        planetarium->PlotMethod2(celestial_trajectory,
                                 /*first_time=*/plugin->CurrentTime(),
                                 /*last_time=*/final_time,
                                 /*now=*/plugin->CurrentTime(),
                                 /*reverse=*/false);
    return m.Return(new TypedIterator<RP2Lines<Length, Camera>>(rp2_lines));
  }
}

}  // namespace interface
}  // namespace principia
