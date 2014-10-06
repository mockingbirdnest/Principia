
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/10/05-19:16:12
// Benchmark                                  Time(ns)    CPU(ns) Iterations
// -------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly            24203858590 23649751600          1                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly            24101375556 24024154000          1                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly            24065379835 23696551900          1                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly            24012373186 23712152000          1                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly            24073379016 23712152000          1                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean       24091273237 23758952300          1                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev       63235130  134559341          0                                 1.0002759262592109e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies        50587372995 49951520200          1                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies        51049073196 50403923100          1                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies        50057970024 49343116300          1                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies        50023569084 49046714400          1                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies        51184083438 50653524700          1                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean   50580413747 49879759740          1                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev  483141931  610009791          0                                 1.0002759263053576e+00 ua  // NOLINT(whitespace/line_length)
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

BENCHMARK(BM_SolarSystemMajorBodiesOnly);
BENCHMARK(BM_SolarSystemMinorAndMajorBodies);

}  // namespace benchmarks
}  // namespace principia
