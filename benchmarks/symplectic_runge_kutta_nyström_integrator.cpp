
// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>
#include <functional>
#include <type_traits>
#include <vector>

#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"
#include "testing_utilities/integration.hpp"

namespace principia {

using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::Velocity;
using integrators::CompositionMethod;
using integrators::IntegrationProblem;
using integrators::SpecialSecondOrderDifferentialEquation;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::si::Metre;
using testing_utilities::ComputeHarmonicOscillatorAcceleration;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace benchmarks {

namespace {

using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, true>;

using ODE = SpecialSecondOrderDifferentialEquation<Position<World>>;

// TODO(egg): use the one from testing_utilities/integration again when everyone
// uses |Instant|s.
// Increments |*evaluations| if |evaluations| is not null.
void ComputeHarmonicOscillatorAcceleration(
    Instant const& t,
    std::vector<Position<World>> const& q,
    std::vector<Vector<Acceleration, World>>* const result,
    int* evaluations) {
  (*result)[0] =
      (World::origin - q[0]) * (SIUnit<Stiffness>() / SIUnit<Mass>());
  if (evaluations != nullptr) {
    ++*evaluations;
  }
}

}  // namespace

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError(
    not_null<benchmark::State*> const state,
    not_null<Length*> const q_error,
    not_null<Speed*> const v_error,
    Integrator const& integrator) {
  Displacement<World> const q_initial({1 * Metre, 0 * Metre, 0 * Metre});
  Velocity<World> const v_initial;
  Instant const t_initial;
#ifdef _DEBUG
  Instant const t_final = t_initial + 100 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 3.0E-4 * Second;

  int evaluations = 0;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::SystemState const initial_state = {{World::origin + q_initial},
                                          {v_initial},
                                          t_initial};
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  problem.append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  integrator.Solve(problem, step);

  state->PauseTiming();
  *q_error = Length();
  *v_error = Speed();
  for (std::size_t i = 0; i < solution.size(); ++i) {
  auto x = (solution[i].positions[0].value - World::origin);
    *q_error = std::max(*q_error,
                        ((solution[i].positions[0].value - World::origin) -
                            q_initial *
                            Cos((solution[i].time.value - t_initial) *
                                (Radian / Second))).Norm());
    *v_error = std::max(*v_error,
                        (solution[i].velocities[0].value +
                            (q_initial / Second) *
                            Sin((solution[i].time.value - t_initial) *
                                (Radian / Second))).Norm());
  }
  state->ResumeTiming();
}

template<typename Integrator, Integrator const& (*integrator)()>
void BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator(
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

BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::McLachlanAtela1992Order4Optimal<Position<World>>()),
    &integrators::McLachlanAtela1992Order4Optimal<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::McLachlan1995SB3A4<Position<World>>()),
    &integrators::McLachlan1995SB3A4<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::McLachlan1995SB3A5<Position<World>>()),
    &integrators::McLachlan1995SB3A5<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::BlanesMoan2002SRKN6B<Position<World>>()),
    &integrators::BlanesMoan2002SRKN6B<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::McLachlanAtela1992Order5Optimal<Position<World>>()),
    &integrators::McLachlanAtela1992Order5Optimal<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::OkunborSkeel1994Order6Method13<Position<World>>()),
    &integrators::OkunborSkeel1994Order6Method13<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::BlanesMoan2002SRKN11B<Position<World>>()),
    &integrators::BlanesMoan2002SRKN11B<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator,
    decltype(integrators::BlanesMoan2002SRKN14A<Position<World>>()),
    &integrators::BlanesMoan2002SRKN14A<Position<World>>);

}  // namespace benchmarks
}  // namespace principia
