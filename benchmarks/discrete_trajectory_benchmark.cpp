// .\Release\x64\benchmarks.exe --benchmark_filter=BM_DiscreteTrajectory --benchmark_repetitions=5  // NOLINT(whitespace/line_length)

#include <optional>
#include <vector>

#include "astronomy/time_scales.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "benchmark/benchmark.h"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "glog/logging.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using namespace principia::astronomy::_time_scales;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_segment;
using namespace principia::physics::_discrete_trajectory_types;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_discrete_trajectory_factories;
using namespace principia::testing_utilities::_solar_system_factory;

namespace {

Ephemeris<Barycentric>::FixedStepParameters EphemerisParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymmetricLinearMultistepIntegrator<
          QuinlanTremaine1990Order12,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/10 * Minute);
}

Ephemeris<Barycentric>::FixedStepParameters HistoryParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymmetricLinearMultistepIntegrator<
          Quinlan1999Order8A,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/10 * Second);
}

// Constructs a trajectory by assigning the points in `timeline` to segments
// defined by `splits`, which must be doubles in [0, 1].
DiscreteTrajectory<World> MakeTrajectory(Timeline<World> const& timeline,
                                         std::vector<double> const& splits) {
  DiscreteTrajectory<World> trajectory;
  Instant const t_min = timeline.begin()->time;
  Instant const t_max = timeline.rbegin()->time;
  auto it_split = splits.begin();
  std::optional<Instant> t_split;
  for (auto const& [t, degrees_of_freedom] : timeline) {
    // The computation of `t_split` is complicated enough that we want to do it
    // at a single place.
    if (!t_split.has_value()) {
      if (it_split == splits.end()) {
        t_split = InfiniteFuture;
      } else {
        double const split = *it_split;
        CHECK_LE(0, split);
        CHECK_LE(split, 1);
        t_split = Barycentre({t_min, t_max}, {1 - split, split});
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

void BM_DiscreteTrajectoryDownsampling(benchmark::State& state) {
  SolarSystem<Barycentric> const solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt",
      /*ignore_frame=*/true);
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      EphemerisParameters());
  auto const earth = solar_system.massive_body(
      *ephemeris, SolarSystemFactory::name(SolarSystemFactory::Earth));

  // Two-line elements for GOES-8:
  // 1 23051U 94022A   00004.06628221 -.00000243  00000-0  00000-0 0  9630
  // 2 23051   0.4232  97.7420 0004776 192.8349 121.5613  1.00264613 28364
  constexpr Instant goes_8_epoch = "JD2451548.56628221"_UT1;
  KeplerianElements<Barycentric> goes_8_elements;
  goes_8_elements.inclination = 0.4232 * Degree;
  goes_8_elements.longitude_of_ascending_node = 97.7420 * Degree;
  goes_8_elements.eccentricity = 0.0004776;
  goes_8_elements.argument_of_periapsis = 192.8349 * Degree;
  goes_8_elements.mean_anomaly = 121.5613 * Degree;
  goes_8_elements.mean_motion = 1.00264613 * (2 * π * Radian / Day);
  CHECK_OK(ephemeris->Prolong(goes_8_epoch));
  KeplerOrbit<Barycentric> const goes_8_orbit(
      *earth, MasslessBody{}, goes_8_elements, goes_8_epoch);

  // First build a realistic trajectory without downsampling.
  DiscreteTrajectory<Barycentric> goes_8_trajectory;
  CHECK_OK(goes_8_trajectory.Append(
      goes_8_epoch,
      ephemeris->trajectory(earth)->EvaluateDegreesOfFreedom(goes_8_epoch) +
          goes_8_orbit.StateVectors(goes_8_epoch)));
  auto goes_8_instance =
      ephemeris->NewInstance({&goes_8_trajectory},
                             Ephemeris<Barycentric>::NoIntrinsicAccelerations,
                             HistoryParameters());
  CHECK_OK(
      ephemeris->FlowWithFixedStep(goes_8_epoch + 100 * Day, *goes_8_instance));

  // Now append the same points to a trajectory with downsampling.
  DiscreteTrajectory<Barycentric> downsampled_trajectory;
  for (auto _ : state) {
    downsampled_trajectory.clear();
    downsampled_trajectory.segments().begin()->SetDownsampling(
        DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters{
            .max_dense_intervals = 10'000, .tolerance = 10 * Metre});
    for (auto const& [t, degrees_of_freedom] : goes_8_trajectory) {
      CHECK_OK(downsampled_trajectory.Append(t, degrees_of_freedom));
    }
  }
  state.SetLabel((std::stringstream()
                  << goes_8_trajectory.size() << " points before downsampling, "
                  << downsampled_trajectory.size() << " after")
                     .str());
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
BENCHMARK(BM_DiscreteTrajectoryDownsampling);

}  // namespace physics
}  // namespace principia
