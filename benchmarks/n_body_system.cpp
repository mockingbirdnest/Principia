
// .\Release\benchmarks.exe  --benchmark_filter=Solar --benchmark_min_time=120 --benchmark_repetitions=5                // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/11-21:43:28
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        29550253368 29499789100          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        29623949266 29562189500          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        29542950130 29484189000          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        29571937299 29515389200          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        29341910577 29328188000          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   29526200128 29477948960          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev   96403882   79300082          0                                 1.0002759262591376e+000 ua    // NOLINT(whitespace/line_length)
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
