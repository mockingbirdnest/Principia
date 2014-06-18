
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/15-23:58:26
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        23449882246 23400150000          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        24191402412 24164554900          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        22669252611 22651345200          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        22655250050 22604544900          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        23051288582 22963347200          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   23203415180 23156788440          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  573728498  578391459          0                                 1.0002759262590839e+000 ua    // NOLINT(whitespace/line_length)
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
