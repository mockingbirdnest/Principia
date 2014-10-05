
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/10/05-13:01:28
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        48122102714 47954707400          1                                 1.0002759263053012e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        48154779410 48126308500          1                                 1.0002759263053012e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        52067173219 52041933600          1                                 1.0002759263053012e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        48369799829 48188708900          1                                 1.0002759263053012e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        48255792594 48079508200          1                                 1.0002759263053012e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   48993929553 48878233320          1                                 1.0002759263053012e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev 1539055926 1583710596          0                                 1.0002759263053012e+00 ua    // NOLINT(whitespace/line_length)
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

void BM_SolarSystem(benchmark::State& state) {  // NOLINT(runtime/references)
  std::vector<quantities::Momentum> output;
  while (state.KeepRunning()) {
    state.PauseTiming();
    std::unique_ptr<SolarSystem> solar_system = SolarSystem::AtСпутник1Launch();
    state.ResumeTiming();
    SimulateSolarSystem(solar_system.get());
    state.PauseTiming();
    state.SetLabel(
        DebugString(
            (solar_system->trajectories()[0]->last_position() -
             solar_system->trajectories()[5]->last_position()).Norm() /
                AstronomicalUnit) + " ua");
    state.ResumeTiming();
  }
}
BENCHMARK(BM_SolarSystem);

}  // namespace benchmarks
}  // namespace principia
