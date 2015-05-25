
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <vector>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using quantities::Abs;
using quantities::Length;
using si::Centi;
using si::Metre;
using si::Milli;
using si::Second;
using testing_utilities::AbsoluteError;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;

namespace integrators {

using ODE = SpecialSecondOrderDifferentialEquation<Length>;

namespace {

// TODO(egg): use the one from testing_utilities/integration again when everyone
// uses |Instant|s.
void ComputeHarmonicOscillatorAcceleration(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result,
    not_null<int*> evaluations) {
  (*result)[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
  ++*evaluations;
}

double HarmonicOscillatorToleranceRatio(
    Time const& h,
    ODE::SystemStateError const& error,
    Length const& q_tolerance,
    Speed const& v_tolerance,
    std::function<void(bool tolerable)> callback) {
  double const r = std::min(q_tolerance / Abs(error.position_error[0]),
                            v_tolerance / Abs(error.velocity_error[0]));
  callback(r > 1.0);
  return r;
}

}  // namespace

class EmbeddedExplicitRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest,
       HarmonicOscillatorBackAndForth) {
  AdaptiveSizeIntegrator<ODE> const& integrator =
      DormandElMikkawyPrince1986RKN434FM<Length>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;
  int const steps_forward = 132;
  // We integrate backward with double the tolerance.
  int const steps_backward = 112;

  int evaluations = 0;
  int initial_rejections = 0;
  int subsequent_rejections = 0;
  bool first_step = true;
  auto const step_size_callback = [&initial_rejections, &subsequent_rejections,
                                   &first_step](bool tolerable) {
    if (!tolerable) {
      if (first_step) {
        ++initial_rejections;
      } else {
        ++subsequent_rejections;
      }
    } else if (first_step) {
      first_step = false;
    }
  };

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::SystemState const initial_state = {{x_initial}, {v_initial}, t_initial};
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  problem.append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };
  AdaptiveStepSize<ODE> adaptive_step_size;
  adaptive_step_size.first_time_step = t_final - t_initial;
  adaptive_step_size.safety_factor = 0.9;
  adaptive_step_size.tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, length_tolerance, speed_tolerance, step_size_callback);

  integrator.Solve(problem, adaptive_step_size);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(3E-4 * Metre), Le(4E-4 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(2E-3 * Metre / Second), Le(3E-3 * Metre / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  EXPECT_EQ(steps_forward, solution.size());
  EXPECT_EQ(
      (1 + initial_rejections) * 4 +
          (steps_forward - 1 + subsequent_rejections) * 3,
      evaluations);

  evaluations = 0;
  subsequent_rejections = 0;
  initial_rejections = 0;
  first_step = true;
  problem.initial_state = &solution.back();
  problem.t_final = t_initial;
  adaptive_step_size.first_time_step = t_initial - t_final;
  adaptive_step_size.tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, 2 * length_tolerance, 2 * speed_tolerance,
                step_size_callback);

  integrator.Solve(problem, adaptive_step_size);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(1E-3 * Metre), Le(2E-3 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(2E-3 * Metre / Second), Le(3E-3 * Metre / Second)));
  EXPECT_EQ(t_initial, solution.back().time.value);
  EXPECT_EQ(steps_backward, solution.size() - steps_forward);
  EXPECT_EQ(
      (1 + initial_rejections) * 4 +
          (steps_backward - 1 + subsequent_rejections) * 3,
      evaluations);
}

}  // namespace integrators
}  // namespace principia
