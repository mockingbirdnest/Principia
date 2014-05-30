
// .\Release\benchmarks_tests.exe --benchmark_repetitions=5 --benchmark_min_time=300
// Benchmarking on 1 X 3310 MHz CPU
// 2014/05/30-17:19:56
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_SolveHarmonicOscillator        5503073850 5488398818         11                                 1.37019e-013, 1.37057e-013
// BM_SolveHarmonicOscillator        5066959920 5054432400         12                                 1.37019e-013, 1.37057e-013
// BM_SolveHarmonicOscillator        5647669798 5640145245         11                                 1.37019e-013, 1.37057e-013
// BM_SolveHarmonicOscillator        5107097166 5098832685         13                                 1.37019e-013, 1.37057e-013
// BM_SolveHarmonicOscillator        4840442968 4828830954         13                                 1.37019e-013, 1.37057e-013
// BM_SolveHarmonicOscillator_mean   5212995348 5202113347         60                                 1.37019e-013, 1.37057e-013
// BM_SolveHarmonicOscillator_stddev  294675070  295069109         60                                 1.37019e-013, 1.37057e-013

//C:\Users\phl\Projects\GitHub\Principia [Tuning]> .\Release\benchmarks_tests.exe --benchmark_repetitions=5 --benchmark_min_time=300
//Benchmarking on 1 X 3310 MHz CPU
//2014/05/30-18:06:19
//Benchmark                           Time(ns)    CPU(ns) Iterations
//------------------------------------------------------------------
//BM_SolveHarmonicOscillator        2320474628 2319614869         26                                 1.37019e-013, 1.37057e-013
//BM_SolveHarmonicOscillator        2241096931 2237169896         27                                 1.37019e-013, 1.37057e-013
//BM_SolveHarmonicOscillator        2146287530 2143130979         29                                 1.37019e-013, 1.37057e-013
//BM_SolveHarmonicOscillator        2268943885 2267214533         27                                 1.37019e-013, 1.37057e-013
//BM_SolveHarmonicOscillator        2448445449 2442975660         25                                 1.37019e-013, 1.37057e-013
//BM_SolveHarmonicOscillator_mean   2280275500 2277265344        134                                 1.37019e-013, 1.37057e-013
//BM_SolveHarmonicOscillator_stddev   99047135   98535204        134                                 1.37019e-013, 1.37057e-013

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

  parameters.q0 = {1.0};
  parameters.p0 = {0.0};
  parameters.t0 = 0.0;
#ifdef _DEBUG
  parameters.tmax = 100.0;
#else
  parameters.tmax = 1000.0;
#endif
  parameters.Δt = 1.0E-4;
  parameters.coefficients = integrator.Order5Optimal();
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
