
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/10/08-21:52:53
// Benchmark                                     Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly               22719357706 22604544900          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22830249525 22776146000          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22695240236 22620145000          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22681241251 22620145000          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22631234623 22526544400          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean          22711464668 22629505060          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev            66008620   81120520          0                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           47008883453 46878300500          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           52003164268 51839132300          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           51617127157 51449129800          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           47164680220 47034301500          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46987666822 46956301000          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean      48956304384 48831433020          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev      2334147907 2300396480          0                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        59040446020 58968378000          1                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        58897055841 58859177300          1                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        58945859648 58905977600          1                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        58946855522 58921577700          1                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        63729339338 63710808400          1                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_mean   59911911274 59873183800          1                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_stddev   1909278849 1919129349          0                                 1.0002759263000629e+00 ua  // NOLINT(whitespace/line_length)
#include <memory>
#include <vector>

#include "benchmarks/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::testing_utilities::ICRFJ2000Ecliptic;
using principia::physics::NBodySystem;
using principia::quantities::DebugString;
using principia::si::AstronomicalUnit;

namespace principia {
namespace benchmarks {

namespace {

void SolarSystemBenchmark(SolarSystem::Accuracy const accuracy,
                          benchmark::State* state) {
  std::vector<quantities::Momentum> output;
  while (state->KeepRunning()) {
    state->PauseTiming();
    std::unique_ptr<SolarSystem> solar_system = SolarSystem::AtСпутник1Launch(
        accuracy);
    state->ResumeTiming();
    SimulateSolarSystem(solar_system.get());
    state->PauseTiming();
    state->SetLabel(
        DebugString(
            (solar_system->trajectories()[
                 SolarSystem::kSun]->last_position() -
             solar_system->trajectories()[
                 SolarSystem::kEarth]->last_position()).Norm() /
                 AstronomicalUnit) + " ua");
    state->ResumeTiming();
  }
}

}  // namespace

void BM_SolarSystemMajorBodiesOnly(
    benchmark::State& state) {  // NOLINT(runtime/references)
  SolarSystemBenchmark(SolarSystem::Accuracy::kMajorBodiesOnly,
                       &state);
}

void BM_SolarSystemMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  SolarSystemBenchmark(SolarSystem::Accuracy::kMinorAndMajorBodies,
                       &state);
}

void BM_SolarSystemAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  SolarSystemBenchmark(SolarSystem::Accuracy::kAllBodiesAndOblateness,
                       &state);
}

BENCHMARK(BM_SolarSystemMajorBodiesOnly);
BENCHMARK(BM_SolarSystemMinorAndMajorBodies);
BENCHMARK(BM_SolarSystemAllBodiesAndOblateness);

}  // namespace benchmarks
}  // namespace principia
