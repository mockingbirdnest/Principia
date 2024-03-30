// .\Release\x64\benchmarks.exe --benchmark_repetitions=5 --benchmark_filter=EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>
#include <functional>
#include <vector>

#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/integration.hpp"

namespace principia {
namespace integrators {

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_integration;

namespace {

using World = Frame<struct WorldTag, Inertial>;

using ODE1D = SpecialSecondOrderDifferentialEquation<Length>;
using ODE3D = SpecialSecondOrderDifferentialEquation<Position<World>>;

template<typename ODE>
double HarmonicOscillatorToleranceRatio1D(
    Time const& h,
    typename ODE::State const& /*state*/,
    typename ODE::State::Error const& error,
    Length const& q_tolerance,
    Speed const& v_tolerance) {
  return std::min(q_tolerance / Abs(error.position_error[0]),
                  v_tolerance / Abs(error.velocity_error[0]));
}

template<typename ODE>
double HarmonicOscillatorToleranceRatio3D(
    Time const& h,
    typename ODE::State const& /*state*/,
    typename ODE::State::Error const& error,
    Length const& q_tolerance,
    Speed const& v_tolerance) {
  return std::min(q_tolerance / (error.position_error[0]).Norm(),
                  v_tolerance / (error.velocity_error[0]).Norm());
}

}  // namespace

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError1D(benchmark::State& state,
                                              Length& q_error,
                                              Speed& v_error,
                                              Integrator const& integrator) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial;
  Instant const t_initial;
  Instant const t_final = t_initial + 1000 * Second;
  Length const length_tolerance = 1e-6 * Metre;
  Speed const speed_tolerance = 1e-6 * Metre / Second;

  std::vector<ODE1D::State> solution;
  ODE1D harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE1D> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE1D::State const& state) {
    solution.push_back(state);
  };

  typename Integrator::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio1D<ODE1D>,
                _1, _2, _3, length_tolerance, speed_tolerance);

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  CHECK_OK(instance->Solve(t_final));

  state.PauseTiming();
  q_error = Length();
  v_error = Speed();
  for (std::size_t i = 0; i < solution.size(); ++i) {
    q_error = std::max(q_error,
                       Abs(solution[i].positions[0].value  -
                           Metre *
                           Cos((solution[i].time.value - t_initial) *
                               (Radian / Second))));
    v_error = std::max(v_error,
                       Abs(solution[i].velocities[0].value +
                           (Metre / Second) *
                           Sin((solution[i].time.value - t_initial) *
                               (Radian / Second))));
  }
  state.ResumeTiming();
}

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError3D(
    benchmark::State& state,
    Length& q_error,
    Speed& v_error,
    Integrator const& integrator) {
  Displacement<World> const q_initial({1 * Metre, 0 * Metre, 0 * Metre});
  Velocity<World> const v_initial;
  Instant const t_initial;
  Instant const t_final = t_initial + 1000 * Second;
  Length const length_tolerance = 1e-6 * Metre;
  Speed const speed_tolerance = 1e-6 * Metre / Second;

  std::vector<ODE3D::State> solution;
  ODE3D harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration3D<World>,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE3D> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {World::origin + q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE3D::State const& state) {
    solution.push_back(state);
  };

  typename Integrator::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio3D<ODE3D>,
                _1, _2, _3, length_tolerance, speed_tolerance);

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  CHECK_OK(instance->Solve(t_final));

  state.PauseTiming();
  q_error = Length();
  v_error = Speed();
  for (std::size_t i = 0; i < solution.size(); ++i) {
    q_error = std::max(q_error,
                       ((solution[i].positions[0].value - World::origin) -
                           q_initial *
                           Cos((solution[i].time.value - t_initial) *
                               (Radian / Second))).Norm());
    v_error = std::max(v_error,
                       (solution[i].velocities[0].value +
                           (q_initial / Second) *
                           Sin((solution[i].time.value - t_initial) *
                               (Radian / Second))).Norm());
  }
  state.ResumeTiming();
}

template<typename Method, typename ODE>
void BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  for (auto _ : state) {
    SolveHarmonicOscillatorAndComputeError1D(
        state,
        q_error,
        v_error,
        EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE>());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

template<typename Method, typename ODE>
void BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  for (auto _ : state) {
    SolveHarmonicOscillatorAndComputeError3D(
        state,
        q_error,
        v_error,
        EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, ODE>());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

// Keep each argument on a single line below, lest it breaks benchmark parsing.

BENCHMARK_TEMPLATE2(
    BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::DormandالمكاوىPrince1986RKN434FM, ODE1D)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE2(
    BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::DormandالمكاوىPrince1986RKN434FM, ODE3D)
    ->Unit(benchmark::kMillisecond);

}  // namespace integrators
}  // namespace principia
