
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/10/09-22:16:45
// Benchmark                                     Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly               22727676130 22401743600          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               23477306582 23228548900          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22756247021 22682545400          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean          22987076578 22770945967          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev          346841142   343280979          0                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46018670298 45770693400          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46508613325 46347897100          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           44590424753 44444684900          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean      45705902792 45521091800          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev      813727835   796776794          0                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56108044602 55988758900          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56511614538 56269560700          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56057569719 55848358000          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_mean   56225742953 56035559200          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_stddev   203189322   175110688          0                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
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
