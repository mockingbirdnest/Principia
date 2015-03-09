
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2015/03/09-23:32:30
// Benchmark                                     Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly               22742077566 22588944800          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               23056272484 22807346200          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22857255198 22620145000          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean          22885201749 22672145333          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev          129782816    96446251          0                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46631853558 46410297500          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           46360598541 46176296000          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           48228785969 47892307000          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean      47073746023 46826300167          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev      824209835   759810109          0                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56597427584 56472362000          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56928653933 56550362500          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        56740186430 56550362500          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_mean   56755422649 56524362333          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_stddev   135651097    36769788          0                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)

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
