
// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=VisibleSegments  // NOLINT(whitespace/line_length)

#include "geometry/perspective.hpp"

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

using quantities::Length;
using quantities::si::Metre;

namespace {
using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST1, false>;
using Camera = Frame<serialization::Frame::TestTag,
                     serialization::Frame::TEST2, false>;
}  // namespace

void RandomSegmentsBenchmark(
    std::uniform_real_distribution<> const& x_distribution,
    std::uniform_real_distribution<> const& y_distribution,
    std::uniform_real_distribution<> const& z_distribution,
    benchmark::State& state) {
  // The camera is on the x-axis and looks towards the positive x.
  Point<Displacement<World>> const camera_origin(
      World::origin +
      Displacement<World>({-10 * Metre, 0 * Metre, 0 * Metre}));
  AffineMap<World, Camera, Length, OrthogonalMap> const world_to_camera_affine(
      camera_origin,
      Camera::origin,
      OrthogonalMap<World, Camera>::Identity());
  Perspective<World, Camera, Length, OrthogonalMap> const perspective(
      world_to_camera_affine,
      /*focal=*/1 * Metre);

  // The sphere is at the origin and has unit radius.
  Sphere<Length, World> const sphere(World::origin,
                                     /*radius=*/1 * Metre);

  int const count = state.range_x();
  std::mt19937_64 random(42);
  std::vector<Segment<Displacement<World>>> segments;
  for (int i = 0; i < count; ++i) {
    segments.emplace_back(
        Point<Displacement<World>>(
            World::origin +
            Displacement<World>({x_distribution(random) * Metre,
                                 y_distribution(random) * Metre,
                                 z_distribution(random) * Metre})),
        Point<Displacement<World>>(
            World::origin +
            Displacement<World>({x_distribution(random) * Metre,
                                 y_distribution(random) * Metre,
                                 z_distribution(random) * Metre})));
  }

  int visible_segments_count = 0;
  int visible_segments_size = 0;
  while (state.KeepRunning()) {
    for (auto const& segment : segments) {
      auto const visible_segments =
          perspective.VisibleSegments(segment, sphere);
      ++visible_segments_count;
      visible_segments_size += visible_segments.size();
    };
  }

  state.SetLabel("average visible segments: " +
                 std::to_string(static_cast<double>(visible_segments_size) /
                                static_cast<double>(visible_segments_count)));
}

void BM_VisibleSegmentsRandomEverywhere(benchmark::State& state) {
  // Generate random segments in the cube [-10, 10[³.
  std::uniform_real_distribution<> distribution(-10.0, 10.0);
  RandomSegmentsBenchmark(distribution, distribution, distribution, state);
}

void BM_VisibleSegmentsRandomNoIntersection(benchmark::State& state) {
  // Generate random segments in the volume [-10, 10[² × [10, 20[.  Note that
  // there is no guarantee that the sphere will not occasionally intersect the
  // plane KAP
  std::uniform_real_distribution<> xy_distribution(-10.0, 10.0);
  std::uniform_real_distribution<> z_distribution(10.0, 20.0);
  RandomSegmentsBenchmark(
      xy_distribution, xy_distribution, z_distribution, state);
}

BENCHMARK(BM_VisibleSegmentsRandomEverywhere)->Arg(1000);
BENCHMARK(BM_VisibleSegmentsRandomNoIntersection)->Arg(1000);

}  // namespace geometry
}  // namespace principia
