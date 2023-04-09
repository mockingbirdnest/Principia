// .\Release\x64\benchmarks.exe --benchmark_filter=DiscreteTrajectory --benchmark_repetitions=5  // NOLINT(whitespace/line_length)

#include "physics/discrete_trajectory.hpp"

#include <vector>

#include "astronomy/epoch.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"

namespace principia {
namespace physics {

using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_discrete_trajectory_types;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_discrete_trajectory_factories;

namespace {

// Constructs a trajectory by assigning the points in |timeline| to segments
// defined by |splits|, which must be doubles in [0, 1].
DiscreteTrajectory<World> MakeTrajectory(Timeline<World> const& timeline,
                                         std::vector<double> const& splits) {
  DiscreteTrajectory<World> trajectory;
  Instant const t_min = timeline.begin()->time;
  Instant const t_max = timeline.rbegin()->time;
  auto it_split = splits.begin();
  std::optional<Instant> t_split;
  for (auto const& [t, degrees_of_freedom] : timeline) {
    // The computation of |t_split| is complicated enough that we want to do it
    // at a single place.
    if (!t_split.has_value()) {
      if (it_split == splits.end()) {
        t_split = InfiniteFuture;
      } else {
        double const split = *it_split;
        CHECK_LE(0, split);
        CHECK_LE(split, 1);
        t_split = Barycentre<Instant, double>({t_min, t_max},
                                              {1 - split, split});
      }
    }
    if (t >= t_split) {
      ++it_split;
      t_split = std::nullopt;
      trajectory.NewSegment();
    }
    CHECK_OK(trajectory.Append(t, degrees_of_freedom));
  }
  return trajectory;
}

// Constructs a trajectory beginning with many empty segments.
DiscreteTrajectory<World> MakeTrajectoryWithEmptySegments(
    int const number_of_empty_segments) {
  DiscreteTrajectory<World> trajectory;
  for (int i = 0; i < number_of_empty_segments; i++) {
    trajectory.NewSegment();
  }
  CHECK_OK(trajectory.Append(
      Instant(), DegreesOfFreedom<World>(World::origin, Velocity<World>())));
  return trajectory;
}

}  // namespace

void BM_DiscreteTrajectoryFront(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.front());
  }
}

void BM_DiscreteTrajectoryFrontEmpty(benchmark::State& state) {
  auto const trajectory = MakeTrajectoryWithEmptySegments(100);

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.front());
  }
}

void BM_DiscreteTrajectoryBack(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.back());
  }
}

void BM_DiscreteTrajectoryBegin(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.begin());
  }
}

void BM_DiscreteTrajectoryBeginEmpty(benchmark::State& state) {
  auto const trajectory = MakeTrajectoryWithEmptySegments(100);

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.begin());
  }
}

void BM_DiscreteTrajectoryEnd(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.end());
  }
}

void BM_DiscreteTrajectoryTMin(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.t_min());
  }
}

void BM_DiscreteTrajectoryTMax(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    benchmark::DoNotOptimize(trajectory.t_max());
  }
}

void BM_DiscreteTrajectorySegmentFront(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.front());
  }
}

void BM_DiscreteTrajectorySegmentFrontEmpty(benchmark::State& state) {
  auto const trajectory = MakeTrajectoryWithEmptySegments(100);
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.front());
  }
}

void BM_DiscreteTrajectorySegmentBack(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.back());
  }
}

void BM_DiscreteTrajectorySegmentBegin(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.begin());
  }
}

void BM_DiscreteTrajectorySegmentBeginEmpty(benchmark::State& state) {
  auto const trajectory = MakeTrajectoryWithEmptySegments(100);
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.begin());
  }
}

void BM_DiscreteTrajectorySegmentEnd(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.end());
  }
}

void BM_DiscreteTrajectorySegmentTMin(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.t_min());
  }
}

void BM_DiscreteTrajectorySegmentTMax(benchmark::State& state) {
  Instant const t0;
  auto const timeline = NewMotionlessTrajectoryTimeline(World::origin,
                                                        /*Δt=*/1 * Second,
                                                        /*t1=*/t0,
                                                        /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});
  auto const& segment = trajectory.segments().front();

  for (auto _ : state) {
    benchmark::DoNotOptimize(segment.t_max());
  }
}

void BM_DiscreteTrajectoryCreateDestroy(benchmark::State& state) {
  Instant const t0;
  int const steps = state.range(0);
  auto const timeline =
      NewMotionlessTrajectoryTimeline(World::origin,
                                      /*Δt=*/1 * Second,
                                      /*t1=*/t0,
                                      /*t2=*/t0 + steps * Second);
  for (auto _ : state) {
    MakeTrajectory(timeline, {0.5, 0.75});
  }
}

void BM_DiscreteTrajectoryIterate(benchmark::State& state) {
  Instant const t0;
  int const steps = state.range(0);
  auto const timeline =
      NewMotionlessTrajectoryTimeline(World::origin,
                                      /*Δt=*/1 * Second,
                                      /*t1=*/t0,
                                      /*t2=*/t0 + steps * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  auto const begin = trajectory.begin();
  auto const end = trajectory.end();
  for (auto _ : state) {
    for (auto it = begin; it != end; ++it) {
    }
  }
}

void BM_DiscreteTrajectoryReverseIterate(benchmark::State& state) {
  Instant const t0;
  int const steps = state.range(0);
  auto const timeline =
      NewMotionlessTrajectoryTimeline(World::origin,
                                      /*Δt=*/1 * Second,
                                      /*t1=*/t0,
                                      /*t2=*/t0 + steps * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  auto const begin = trajectory.begin();
  auto const end = trajectory.end();
  for (auto _ : state) {
    for (auto it = end; it != begin; --it) {
    }
  }
}

void BM_DiscreteTrajectoryFind(benchmark::State& state) {
  Instant const t0;
  int const steps = state.range(0);
  auto const timeline =
      NewMotionlessTrajectoryTimeline(World::origin,
                                      /*Δt=*/1 * Second,
                                      /*t1=*/t0,
                                      /*t2=*/t0 + steps * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    // These times are in different segments of the trajectory.
    trajectory.find(Instant() + 1.0 / 3.0 * steps * Second);
    trajectory.find(Instant() + 2.0 / 3.0 * steps * Second);
    trajectory.find(Instant() + 5.0 / 6.0 * steps * Second);
  }
}

void BM_DiscreteTrajectoryLowerBound(benchmark::State& state) {
  Instant const t0;
  int const steps = state.range(0);
  auto const timeline =
      NewMotionlessTrajectoryTimeline(World::origin,
                                      /*Δt=*/1 * Second,
                                      /*t1=*/t0,
                                      /*t2=*/t0 + steps * Second);
  auto const trajectory = MakeTrajectory(timeline, {0.5, 0.75});

  for (auto _ : state) {
    // These times are in different segments of the trajectory.
    trajectory.lower_bound(Instant() + 1.0 / 3.0 * steps * Second);
    trajectory.lower_bound(Instant() + 2.0 / 3.0 * steps * Second);
    trajectory.lower_bound(Instant() + 5.0 / 6.0 * steps * Second);
  }
}

void BM_DiscreteTrajectoryEvaluateDegreesOfFreedomExact(
    benchmark::State& state) {
  Instant const t0;
  auto const timeline =
      NewCircularTrajectoryTimeline<World>(/*ω=*/1 * Radian / Second,
                                           /*r=*/1 * Metre,
                                           /*Δt=*/1 * Second,
                                           /*t1=*/t0,
                                           /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {});

  Instant const t = Instant() + 1 * Second;
  for (auto _ : state) {
    trajectory.EvaluateDegreesOfFreedom(t);
  }
}

void BM_DiscreteTrajectoryEvaluateDegreesOfFreedomInterpolated(
    benchmark::State& state) {
  Instant const t0;
  auto const timeline =
      NewCircularTrajectoryTimeline<World>(/*ω=*/1 * Radian / Second,
                                           /*r=*/1 * Metre,
                                           /*Δt=*/1 * Second,
                                           /*t1=*/t0,
                                           /*t2=*/t0 + 4 * Second);
  auto const trajectory = MakeTrajectory(timeline, {});

  Instant const t = Instant() + 1.5 * Second;
  for (auto _ : state) {
    trajectory.EvaluateDegreesOfFreedom(t);
  }
}

BENCHMARK(BM_DiscreteTrajectoryFront);
BENCHMARK(BM_DiscreteTrajectoryFrontEmpty);
BENCHMARK(BM_DiscreteTrajectoryBack);
BENCHMARK(BM_DiscreteTrajectoryBegin);
BENCHMARK(BM_DiscreteTrajectoryBeginEmpty);
BENCHMARK(BM_DiscreteTrajectoryEnd);
BENCHMARK(BM_DiscreteTrajectoryTMin);
BENCHMARK(BM_DiscreteTrajectoryTMax);
BENCHMARK(BM_DiscreteTrajectorySegmentFront);
BENCHMARK(BM_DiscreteTrajectorySegmentFrontEmpty);
BENCHMARK(BM_DiscreteTrajectorySegmentBack);
BENCHMARK(BM_DiscreteTrajectorySegmentBegin);
BENCHMARK(BM_DiscreteTrajectorySegmentBeginEmpty);
BENCHMARK(BM_DiscreteTrajectorySegmentEnd);
BENCHMARK(BM_DiscreteTrajectorySegmentTMin);
BENCHMARK(BM_DiscreteTrajectorySegmentTMax);
BENCHMARK(BM_DiscreteTrajectoryCreateDestroy)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryIterate)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryReverseIterate)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryFind)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryLowerBound)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryEvaluateDegreesOfFreedomExact);
BENCHMARK(BM_DiscreteTrajectoryEvaluateDegreesOfFreedomInterpolated);

}  // namespace physics
}  // namespace principia
