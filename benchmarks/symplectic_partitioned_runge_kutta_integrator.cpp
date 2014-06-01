#include "benchmarks/symplectic_partitioned_runge_kutta_integrator.hpp"

// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=300                                                    // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/01-10:56:16
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_SolveHarmonicOscillator        1158462428 1157407419         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1164918290 1163707460         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1161977190 1163107456         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1161909218 1160707440         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1163742772 1162207450         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_mean   1162201980 1161427445        260                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_stddev    2185080    2249814        260                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace benchmarks {

void SolveHarmonicOscillatorAndComputeError(benchmark::State* state,
                                            double* q_error,
                                            double* p_error) {
  SPRKIntegrator::Solution solution;

  SolveHarmonicOscillator(&solution);

  state->PauseTiming();
  *q_error = 0;
  *p_error = 0;
  for (size_t i = 0; i < solution.time.quantities.size(); ++i) {
    *q_error = std::max(*q_error,
                        std::abs(solution.position[0].quantities[i] -
                                 std::cos(solution.time.quantities[i])));
    *p_error = std::max(*p_error,
                        std::abs(solution.momentum[0].quantities[i] +
                                 std::sin(solution.time.quantities[i])));
  }
  state->ResumeTiming();
}

static void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  double q_error;
  double p_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError(&state, &q_error, &p_error);
  }
  std::stringstream ss;
  ss << q_error << ", " << p_error;
  state.SetLabel(ss.str());
}
BENCHMARK(BM_SolveHarmonicOscillator);

}  // namespace benchmarks
}  // namespace principia
