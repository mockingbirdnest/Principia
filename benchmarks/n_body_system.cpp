
// .\Release\benchmarks.exe  --benchmark_filter=Solar --benchmark_min_time=120 --benchmark_repetitions=5                // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/11-21:58:44
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        24890476919 24882159500          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        26308986163 26301768600          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        24715988136 24710558400          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        24851988054 24850959300          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        25200986600 24975760100          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   25193685174 25144241180          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  579799316  584985030          0                                 1.0002759262590839e+000 ua    // NOLINT(whitespace/line_length)
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
