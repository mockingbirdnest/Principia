
// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator                                                 // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/09-00:22:28
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_SolveHarmonicOscillator        1257378650 1260488080          5                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1257609773 1260488080          5                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1260331750 1263608100          5                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1264436937 1263608100          5                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1265400577 1266728120          5                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_mean   1261031538 1262984096          5                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_stddev    3353416    2334809          0                                 1.3701886847350409e-013 m, 1.3705703238997557e-013 m kg s^-1  // NOLINT(whitespace/line_length)

#include "benchmarks/symplectic_partitioned_runge_kutta_integrator.hpp"

#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <algorithm>

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::integrators::SPRKIntegrator;
using principia::quantities::Abs;
using principia::quantities::AngularFrequency;
using principia::quantities::Cos;
using principia::quantities::Length;
using principia::quantities::Momentum;

namespace principia {
namespace benchmarks {

void SolveHarmonicOscillatorAndComputeError(benchmark::State* state,
                                            Length* q_error,
                                            Momentum* p_error) {
  SPRKIntegrator<Length, Momentum>::Solution solution;

  SolveHarmonicOscillator(&solution);

  state->PauseTiming();
  *q_error = Length();
  *p_error = Momentum();
  for (size_t i = 0; i < solution.time.quantities.size(); ++i) {
    *q_error = std::max(*q_error,
                        Abs(solution.position[0].quantities[i] -
                            SIUnit<Length>() *
                            Cos(solution.time.quantities[i] *
                                SIUnit<AngularFrequency>())));
    *p_error = std::max(*p_error,
                        Abs(solution.momentum[0].quantities[i] +
                            SIUnit<Momentum>() *
                            Sin(solution.time.quantities[i] *
                                SIUnit<AngularFrequency>())));
  }
  state->ResumeTiming();
}

static void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  Length   q_error;
  Momentum p_error;
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
