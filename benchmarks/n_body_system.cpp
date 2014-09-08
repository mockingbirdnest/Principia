
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/08/25-22:40:12
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        23378544069 23337749600          1                                 1.0002759262590839e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        23885358072 23852552900          1                                 1.0002759262590839e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        23644338823 23524950800          1                                 1.0002759262590839e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        23363312459 23290949300          1                                 1.0002759262590839e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        23146287656 23103748100          1                                 1.0002759262590839e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   23483568216 23421990140          1                                 1.0002759262590839e+00 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  255551158  253586811          0                                 1.0002759262590839e+00 ua    // NOLINT(whitespace/line_length)
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
    std::unique_ptr<SolarSystem> solar_system = SolarSystem::AtСпутникLaunch();
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
