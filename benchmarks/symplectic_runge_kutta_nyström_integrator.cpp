// .\Release\x64\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=5 --benchmark_filter=SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>
#include <functional>
#include <type_traits>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/methods.hpp"
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
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Mass;
using quantities::Sin;
using quantities::Speed;
using quantities::Stiffness;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::ComputeHarmonicOscillatorAcceleration1D;
using testing_utilities::ComputeHarmonicOscillatorAcceleration3D;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace integrators {

namespace {

using World = Frame<struct WorldTag, Inertial>;
using ODE1D = SpecialSecondOrderDifferentialEquation<Length>;
using ODE3D = SpecialSecondOrderDifferentialEquation<Position<World>>;


}  // namespace

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError1D(benchmark::State& state,
                                              Length& q_error,
                                              Speed& v_error,
                                              Integrator const& integrator) {
  state.PauseTiming();
  Length const q_initial = 1 * Metre;
  Speed const v_initial;
  Instant const t_initial;
#ifdef _DEBUG
  Instant const t_final = t_initial + 100 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 3.0e-4 * Second;

  std::vector<ODE1D::State> solution;
  solution.reserve(static_cast<int>((t_final - t_initial) / step));
  ODE1D harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE1D> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE1D::State const& state) {
    solution.emplace_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);

  state.ResumeTiming();
  CHECK_OK(instance->Solve(t_final));
  state.PauseTiming();

  q_error = Length();
  v_error = Speed();
  for (auto const& state : solution) {
    q_error = std::max(q_error,
                       Abs(state.positions[0].value  -
                           Metre *
                           Cos((state.time.value - t_initial) *
                               (Radian / Second))));
    v_error = std::max(v_error,
                       Abs(state.velocities[0].value +
                           (Metre / Second) *
                           Sin((state.time.value - t_initial) *
                               (Radian / Second))));
  }
  state.ResumeTiming();
}

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError3D(benchmark::State& state,
                                              Length& q_error,
                                              Speed& v_error,
                                              Integrator const& integrator) {
  state.PauseTiming();
  Displacement<World> const q_initial({1 * Metre, 0 * Metre, 0 * Metre});
  Velocity<World> const v_initial;
  Instant const t_initial;
#ifdef _DEBUG
  Instant const t_final = t_initial + 100 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 3.0e-4 * Second;

  std::vector<ODE3D::State> solution;
  solution.reserve(static_cast<int>((t_final - t_initial) / step));
  ODE3D harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration3D<World>,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE3D> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {World::origin + q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE3D::State const& state) {
    solution.emplace_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);

  state.ResumeTiming();
  CHECK_OK(instance->Solve(t_final));
  state.PauseTiming();

  q_error = Length();
  v_error = Speed();
  for (auto const& state : solution) {
    q_error = std::max(q_error,
                       ((state.positions[0].value - World::origin) -
                           q_initial *
                           Cos((state.time.value - t_initial) *
                               (Radian / Second))).Norm());
    v_error = std::max(v_error,
                       (state.velocities[0].value +
                           (q_initial / Second) *
                           Sin((state.time.value - t_initial) *
                               (Radian / Second))).Norm());
  }
  state.ResumeTiming();
}

template<typename Method, typename ODE>
void BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  for (auto _ : state) {
    SolveHarmonicOscillatorAndComputeError1D(
        state,
        q_error,
        v_error,
        SymplecticRungeKuttaNyströmIntegrator<Method, ODE>());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

template<typename Method, typename ODE>
void BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  for (auto _ : state) {
    SolveHarmonicOscillatorAndComputeError3D(
        state,
        q_error,
        v_error,
        SymplecticRungeKuttaNyströmIntegrator<Method, ODE>());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::McLachlanAtela1992Order4Optimal, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::McLachlan1995SB3A4, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::McLachlan1995SB3A5, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::BlanesMoan2002SRKN6B, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::McLachlanAtela1992Order5Optimal, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::OkunborSkeel1994Order6Method13, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::BlanesMoan2002SRKN11B, ODE1D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::BlanesMoan2002SRKN14A, ODE1D)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::McLachlanAtela1992Order4Optimal, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::McLachlan1995SB3A4, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::McLachlan1995SB3A5, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::BlanesMoan2002SRKN6B, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::McLachlanAtela1992Order5Optimal, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::OkunborSkeel1994Order6Method13, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::BlanesMoan2002SRKN11B, ODE3D)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::BlanesMoan2002SRKN14A, ODE3D)
    ->Unit(benchmark::kMillisecond);

}  // namespace integrators
}  // namespace principia
