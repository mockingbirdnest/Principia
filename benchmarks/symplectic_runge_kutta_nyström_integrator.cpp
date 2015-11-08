
// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>
#include <type_traits>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/physics.pb.h"
#include "testing_utilities/integration.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using geometry::Frame;
using integrators::CompositionMethod;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::SIUnit;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;

namespace benchmarks {

template<typename Position, int order_, bool time_reversible_, int evaluations_,
         CompositionMethod composition_>
void SolveHarmonicOscillatorAndComputeError(
    not_null<benchmark::State*> const state,
    not_null<Length*> const q_error,
    not_null<Speed*> const v_error,
    SymplecticRungeKuttaNyströmIntegrator<
        Position, order_, time_reversible_, evaluations_, composition_> const&
            integrator) {
  typename SymplecticRungeKuttaNyströmIntegrator<
      Position, order_, time_reversible_, evaluations_, composition_>::
          Solution<Length, Speed> solution;
  typename SymplecticRungeKuttaNyströmIntegrator<
      Position, order_, time_reversible_, evaluations_, composition_>::
          Parameters<Length, Speed> parameters;

  parameters.initial.positions.emplace_back(SIUnit<Length>());
  parameters.initial.momenta.emplace_back(Speed());
  parameters.initial.time = Time();
#ifdef _DEBUG
  parameters.tmax = 100.0 * SIUnit<Time>();
#else
  parameters.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters.Δt = 1.0E-4 * SIUnit<Time>();
  parameters.sampling_period = 1;
  integrator.SolveTrivialKineticEnergyIncrement<Length>(
      &ComputeHarmonicOscillatorAcceleration,
      parameters,
      &solution);

  state->PauseTiming();
  *q_error = Length();
  *v_error = Speed();
  for (std::size_t i = 0; i < solution.size(); ++i) {
    *q_error = std::max(*q_error,
                        Abs(solution[i].positions[0].value -
                            SIUnit<Length>() *
                            Cos(solution[i].time.value *
                                SIUnit<AngularFrequency>())));
    *v_error = std::max(*v_error,
                        Abs(solution[i].momenta[0].value +
                            SIUnit<Speed>() *
                            Sin(solution[i].time.value *
                                SIUnit<AngularFrequency>())));
  }
  state->ResumeTiming();
}

template<typename Integrator, Integrator const& (*integrator)()>
void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  Length q_error;
  Speed v_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError(&state, &q_error, &v_error,
                                           integrator());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, true>;

BENCHMARK_TEMPLATE2(
    BM_SolveHarmonicOscillator,
    decltype(integrators::McLachlanAtela1992Order4Optimal<Position<World>>()),
    &integrators::McLachlanAtela1992Order4Optimal<Position<World>>);

}  // namespace benchmarks
}  // namespace principia
