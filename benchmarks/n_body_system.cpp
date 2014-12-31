
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// Benchmarking on 1 X 3310 MHz CPU
// 2014/10/14-22:47:11
// Benchmark                                     Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly               23302858330 23306549400          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22130179859 22120941800          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22119180895 22120941800          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean          22517406362 22516144333          1                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev          555416565   558900783          0                                 1.0002759262678920e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46351007201 46332297000          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           45230285860 45177889600          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           45316495873 45318290500          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean      45632596311 45609492367          1                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev      509210956   514304065          0                                 1.0002759263137868e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56332047201 56285160800          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56366599776 56316361000          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56280392863 56253960600          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_mean   56326346613 56285160800          1                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_stddev    35423914    25474857          0                                 1.0002759263001231e+00 ua  // NOLINT(whitespace/line_length)
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
            (solar_system->trajectories()[SolarSystem::kSun]->
                 last().degrees_of_freedom().position() -
             solar_system->trajectories()[SolarSystem::kEarth]->
                 last().degrees_of_freedom().position()).Norm() /
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
