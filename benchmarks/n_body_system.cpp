
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2015/03/10-22:09:40
// Benchmark                                     Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------
// BM_SolarSystemMajorBodiesOnly               22011480785 21824539900          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22581221319 22557744600          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly               22558222271 22495344200          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_mean          22383641458 22292542900          1                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMajorBodiesOnly_stddev          263324786   331907174          0                                 +1.00027592626789200e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           45588232256 45177889600          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           45151473261 44912687900          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies           45157473064 44772287000          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_mean      45299059527 44954288167          1                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemMinorAndMajorBodies_stddev      204490668   168179079          0                                 +1.00027592631378680e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        55712498880 55380355000          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        55391495205 55224354000          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness        55057464815 54974752400          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_mean   55387152967 55193153800          1                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystemAllBodiesAndOblateness_stddev   267434164   167049806          0                                 +1.00027592630012310e+00 ua  // NOLINT(whitespace/line_length)

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
