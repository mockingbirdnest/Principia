
// .\Release\benchmarks.exe --benchmark_filter=Solar --benchmark_min_time=120 --benchmark_repetitions=5                 // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/10-03:59:42
// Benchmark               Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------
// BM_SolarSystem        38965162722 38781848600          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        38893685070 38126644400          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        39232825485 39015850100          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        39306005445 39125050800          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem        38628316847 38485446700          1                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_mean   39005199114 38706968120          5                                 1.0002759262591376e+000 ua  // NOLINT(whitespace/line_length)
// BM_SolarSystem_stddev  244285387  363773460          5                                 1.0002759262591376e+000 ua    // NOLINT(whitespace/line_length)

#include "benchmarks/n_body_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {
namespace benchmarks {

static void BM_SolarSystem(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<quantities::Momentum> output;
  while (state.KeepRunning()) {
    physics::NBodySystem<testing_utilities::ICRFJ2000EclipticFrame>* system =
      SimulateSolarSystem();
    state.PauseTiming();
    state.SetLabel(
        quantities::ToString(
            (system->bodies()[0]->positions().back() -
             system->bodies()[5]->positions().back()).Norm() /
                si::AstronomicalUnit) + " ua");
    state.ResumeTiming();
  }
}
BENCHMARK(BM_SolarSystem);

}  // namespace benchmarks
}  // namespace principia
