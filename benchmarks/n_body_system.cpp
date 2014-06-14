
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/14-19:55:34
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        25027109600 24991360200          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        25547984339 25506163500          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        25297984100 25256561900          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        25230983473 25209761600          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        25151982761 25162961300          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   25251208855 25225361700          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  173631202  166563488          0                                 1.0002759262590839e+000 ua    // NOLINT(whitespace/line_length)
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
