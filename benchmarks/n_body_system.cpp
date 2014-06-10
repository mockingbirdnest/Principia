
// .\Release\benchmarks.exe --benchmark_filter=Solar --benchmark_min_time=120 --benchmark_repetitions=5                 // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/10-05:28:19
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        44898338607 44725486700          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        44167129091 43851881100          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        44233291915 43929881600          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        45446316293 45240290000          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        44918000034 44709886600          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   44732615188 44491485200          5                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  477556136  526812397          5                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)

#include <memory>
#include <vector>

#include "benchmarks/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::testing_utilities::ICRFJ2000EclipticFrame;
using principia::physics::NBodySystem;
using principia::quantities::ToString;
using principia::si::AstronomicalUnit;

namespace principia {
namespace benchmarks {

void BM_SolarSystem(benchmark::State& state) {  // NOLINT(runtime/references)
  std::vector<quantities::Momentum> output;
  while (state.KeepRunning()) {
    std::unique_ptr<NBodySystem<ICRFJ2000EclipticFrame>> const system =
        SimulateSolarSystem();
    state.PauseTiming();
    state.SetLabel(
        ToString(
            (system->bodies()[0]->positions().back() -
             system->bodies()[5]->positions().back()).Norm() /
                AstronomicalUnit) + " ua");
    state.ResumeTiming();
  }
}
BENCHMARK(BM_SolarSystem);

}  // namespace benchmarks
}  // namespace principia
