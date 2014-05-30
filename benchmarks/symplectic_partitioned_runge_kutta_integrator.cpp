
// .\Release\benchmarks_tests.exe --benchmark_repetitions=5 --benchmark_min_time=300                                              // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/05/30-18:44:00
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_SolveHarmonicOscillator        2086500672 2086253373         30                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1835725206 1830411733         33                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        2041165435 2025932987         30                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        2016562322 2011372893         30                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1941378110 1930024872         32                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_mean   1980837962 1973362327        155                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_stddev   88977645   89552186        155                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES
#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include <algorithm>
#include <vector>

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace benchmarks {

namespace {

inline void compute_harmonic_oscillator_force(double const t,
                                              std::vector<double> const& q,
                                              std::vector<double>* result) {
  (*result)[0] = -q[0];
}

inline void compute_harmonic_oscillator_velocity(std::vector<double> const& p,
                                                 std::vector<double>* result) {
  (*result)[0] = p[0];
}

}  // namespace

void SolveHarmonicOscillator(benchmark::State* state,
                             double* q_error,
                             double* p_error) {
  SPRKIntegrator integrator;
  SPRKIntegrator::Parameters parameters;
  SPRKIntegrator::Solution solution;

  integrator.Initialize(integrator.Order5Optimal());

  parameters.q0 = {1.0};
  parameters.p0 = {0.0};
  parameters.t0 = 0.0;
#ifdef _DEBUG
  parameters.tmax = 100.0;
#else
  parameters.tmax = 1000.0;
#endif
  parameters.Δt = 1.0E-4;
  parameters.sampling_period = 1;
  integrator.Solve(&compute_harmonic_oscillator_force,
                   &compute_harmonic_oscillator_velocity,
                   parameters,
                   &solution);

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
    SolveHarmonicOscillator(&state, &q_error, &p_error);
  }
  std::stringstream ss;
  ss << q_error << ", " << p_error;
  state.SetLabel(ss.str());
}
BENCHMARK(BM_SolveHarmonicOscillator);

}  // namespace benchmarks
}  // namespace principia

int main(int argc, const char* argv[]) {
  benchmark::Initialize(&argc, argv);

  benchmark::RunSpecifiedBenchmarks();
}
