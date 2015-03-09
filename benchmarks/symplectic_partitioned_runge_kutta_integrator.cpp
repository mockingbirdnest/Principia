
// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator                                               // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/09/24-00:22:28
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_SolveHarmonicOscillator        3041595356 3010819300          3                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        2969998496 2958818967          3                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        2951685882 2932818800          3                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        2957046406 2948418900          3                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        2960333642 2964019000          3                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_mean   2976131957 2962978993          3                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_stddev   33270202   26186699          0                                 1.3701886847350409e-13 m, 1.3705703238997557e-13 m kg s^-1  // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <vector>

#include "base/not_null.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/numerical_analysis.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using integrators::SPRKIntegrator;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Momentum;
using quantities::SIUnit;
using testing_utilities::ComputeHarmonicOscillatorForce;
using testing_utilities::ComputeHarmonicOscillatorVelocity;

namespace benchmarks {

using Integrator = SPRKIntegrator<Length, Momentum>;

void SolveHarmonicOscillatorAndComputeError(
    not_null<benchmark::State*> const state,
    not_null<Length*> const q_error,
    not_null<Momentum*> const p_error,
    Integrator::Scheme const& (Integrator::*scheme)() const) {
  std::vector<Integrator::SystemState> solution;
  Integrator integrator;
  Integrator::Parameters parameters;

  integrator.Initialize((integrator.*scheme)());

  parameters.initial.positions.emplace_back(SIUnit<Length>());
  parameters.initial.momenta.emplace_back(Momentum());
  parameters.initial.time = Time();
#ifdef _DEBUG
  parameters.tmax = 100.0 * SIUnit<Time>();
#else
  parameters.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters.Δt = 1.0E-4 * SIUnit<Time>();
  parameters.sampling_period = 1;
  integrator.Solve(&ComputeHarmonicOscillatorForce,
                   &ComputeHarmonicOscillatorVelocity,
                   parameters,
                   &solution);

  state->PauseTiming();
  *q_error = Length();
  *p_error = Momentum();
  for (std::size_t i = 0; i < solution.size(); ++i) {
    *q_error = std::max(*q_error,
                        Abs(solution[i].positions[0].value -
                            SIUnit<Length>() *
                            Cos(solution[i].time.value *
                                SIUnit<AngularFrequency>())));
    *p_error = std::max(*p_error,
                        Abs(solution[i].momenta[0].value +
                            SIUnit<Momentum>() *
                            Sin(solution[i].time.value *
                                SIUnit<AngularFrequency>())));
  }
  state->ResumeTiming();
}

template<Integrator::Scheme const& (Integrator::*scheme)() const>
void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  Length   q_error;
  Momentum p_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError(&state, &q_error, &p_error, scheme);
  }
  std::stringstream ss;
  ss << q_error << ", " << p_error;
  state.SetLabel(ss.str());
}

BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Leapfrog);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::PseudoLeapfrog);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order2Optimal);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Ruth1983);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order3Optimal);
BENCHMARK_TEMPLATE(
    BM_SolveHarmonicOscillator,
    &Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta);
BENCHMARK_TEMPLATE(
    BM_SolveHarmonicOscillator,
    &Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order4Optimal);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order5Optimal);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order6A);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order6B);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order6C);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8A);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8B);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8C);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8D);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8E);

}  // namespace benchmarks
}  // namespace principia
