
// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Planetarium  // NOLINT(whitespace/line_length)

#include "ksp_plugin/planetarium.hpp"

#include <random>
#include <vector>

#include "benchmark/benchmark.h"

namespace principia {
namespace geometry {

using base::make_not_null_unique;
using base::not_null;
using geometry::Velocity;
using ksp_plugin::Barycentric;
using ksp_plugin::Planetarium;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using quantities::Angle;
using quantities::Cos;
using quantities::Sin;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

namespace {

not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>>
NewCircularTrajectory(Time const& period, Time const& step, Time const& last) {
  auto discrete_trajectory =
      make_not_null_unique<DiscreteTrajectory<Barycentric>>();
  constexpr Instant t0;
  for (Time t; t <= last; t += step) {
    Angle const ɑ = 2 * π * t * Radian / period;
      DegreesOfFreedom<Barycentric> const degrees_of_freedom(
          Barycentric::origin +
              Displacement<Barycentric>({10 * Metre * Sin(ɑ),
                                         10 * Metre * Cos(ɑ),
                                         0 * Metre}),
          Velocity<Barycentric>({(20 * π / period) * Metre * Cos(ɑ),
                                 -(20 * π / period) * Metre * Sin(ɑ),
                                 0 * Metre / Second}));
      discrete_trajectory->Append(t0 + t, degrees_of_freedom);
  }
  return discrete_trajectory;
}

Planetarium PlanetariumLookingAlongY() {
  // No dark area, human visual acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(
      parameters, perspective_, &ephemeris_, &plotting_frame_);
}

}  // namespace

void BM_PlanetariumPlotMethod2DenseCircleEdgeOn(benchmark::State& state) {
  while (state.KeepRunning()) {
    for (auto const& segment : segments) {
      auto const visible_segments =
          perspective.VisibleSegments(segment, sphere);
      ++visible_segments_count;
      visible_segments_size += visible_segments.size();
    }
  }
}

}  // namespace geometry
}  // namespace principia
