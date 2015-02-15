
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2015/01/10-13:35:55
// Benchmark                                     Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly               22524038053 22479744100          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22509211755 22432943800          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22547224260 22510944300          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean          22526824689 22474544067          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev           15643138    32055158          0                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           47889467931 47471104300          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46522606111 46238696400          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46272590375 46082695400          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean      46894888139 46597498700          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev      710642286   621006783          0                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        57334925151 56503562200          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        57564716077 57189966600          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56853646016 56503562200          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_mean   57251095748 56732363667          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_stddev   296283293   323574137          0                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "benchmarks/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using base::not_null;
using physics::NBodySystem;
using quantities::DebugString;
using si::AstronomicalUnit;
using testing_utilities::ICRFJ2000Ecliptic;

namespace benchmarks {

namespace {

void SolarSystemBenchmark(SolarSystem::Accuracy const accuracy,
                          not_null<benchmark::State*> const state) {
  std::vector<quantities::Momentum> output;
  while (state->KeepRunning()) {
    state->PauseTiming();
    not_null<std::unique_ptr<SolarSystem>> const solar_system =
        SolarSystem::AtСпутник1Launch(accuracy);
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
  SolarSystemBenchmark(SolarSystem::Accuracy::kMajorBodiesOnly, &state);
}

void BM_SolarSystemMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  SolarSystemBenchmark(SolarSystem::Accuracy::kMinorAndMajorBodies, &state);
}

void BM_SolarSystemAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  SolarSystemBenchmark(SolarSystem::Accuracy::kAllBodiesAndOblateness, &state);
}

BENCHMARK(BM_SolarSystemMajorBodiesOnly);
BENCHMARK(BM_SolarSystemMinorAndMajorBodies);
BENCHMARK(BM_SolarSystemAllBodiesAndOblateness);

}  // namespace benchmarks
}  // namespace principia
