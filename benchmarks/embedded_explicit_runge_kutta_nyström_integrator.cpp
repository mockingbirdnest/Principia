
// .\Release\x64\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=5 --benchmark_filter=EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>
#include <functional>
#include <type_traits>
#include <vector>

#include "base/not_null.hpp"
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
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::methods::DormandElMikkawyPrince1986RKN434FM;
using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Mass;
using quantities::Sin;
using quantities::SIUnit;
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

using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, true>;

template<typename ODE>
double HarmonicOscillatorToleranceRatio1D(
    Time const& h,
    typename ODE::SystemStateError const& error,
    Length const& q_tolerance,
    Speed const& v_tolerance) {
  return std::min(q_tolerance / Abs(error.position_error[0]),
                  v_tolerance / Abs(error.velocity_error[0]));
}

template<typename ODE>
double HarmonicOscillatorToleranceRatio3D(
    Time const& h,
    typename ODE::SystemStateError const& error,
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
  using ODE = SpecialSecondOrderDifferentialEquation<Length>;

  Length const q_initial = 1 * Metre;
  Speed const v_initial;
  Instant const t_initial;
  Instant const t_final = t_initial + 1000 * Second;
  Length const length_tolerance = 1e-6 * Metre;
  Speed const speed_tolerance = 1e-6 * Metre / Second;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  Integrator::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio1D<ODE>,
                _1, _2, length_tolerance, speed_tolerance);

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  instance->Solve(t_final);

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
  using ODE = SpecialSecondOrderDifferentialEquation<Position<World>>;

  Displacement<World> const q_initial({1 * Metre, 0 * Metre, 0 * Metre});
  Velocity<World> const v_initial;
  Instant const t_initial;
  Instant const t_final = t_initial + 1000 * Second;
  Time const step = 3.0e-4 * Second;
  Length const length_tolerance = 1e-6 * Metre;
  Speed const speed_tolerance = 1e-6 * Metre / Second;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration3D<World>,
                _1, _2, _3, /*evaluations=*/nullptr);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{World::origin + q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  Integrator::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio3D<ODE>,
                _1, _2, length_tolerance, speed_tolerance);

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  instance->Solve(t_final);

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

template<typename Method, typename Position>
void BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError1D(
        state,
        q_error,
        v_error,
        EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

template<typename Method, typename Position>
void BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError3D(
        state,
        q_error,
        v_error,
        EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

// Keep each argument on a single line below, lest it breaks benchmark parsing.

BENCHMARK_TEMPLATE2(
    BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    methods::DormandElMikkawyPrince1986RKN434FM, Length);

BENCHMARK_TEMPLATE2(
    BM_EmbeddedExplicitRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    methods::DormandElMikkawyPrince1986RKN434FM, Position<World>);

}  // namespace integrators
}  // namespace principia
