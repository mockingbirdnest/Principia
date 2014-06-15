
// .\Release\benchmarks.exe  --benchmark_repetitions=5 --benchmark_filter=Solar
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/15-15:46:17
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        22845606781 22713745600          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        22588242985 22542144500          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        22591244675 22557744600          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        22942278601 22916546900          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        22934276319 22822946300          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   22780329872 22710625580          1                                 1.0002759262590839e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  159272445  146142216          0                                 1.0002759262590839e+000 ua    // NOLINT(whitespace/line_length)
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
