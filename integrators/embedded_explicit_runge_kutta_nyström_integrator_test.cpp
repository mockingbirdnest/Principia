
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <vector>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using quantities::Abs;
using quantities::Length;
using quantities::SpecificImpulse;
using quantities::si::Centi;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Newton;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;
using ::testing::Lt;

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
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
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

  auto outcome = integrator.Solve(problem, adaptive_step_size);
  EXPECT_EQ(TerminationCondition::Done, outcome);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(3e-4 * Metre), Le(4e-4 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(2e-3 * Metre / Second), Le(3e-3 * Metre / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  EXPECT_EQ(steps_forward, solution.size());
  EXPECT_EQ((1 + initial_rejections) * 4 +
                (steps_forward - 1 + subsequent_rejections) * 3,
            evaluations);
  EXPECT_EQ(1, initial_rejections);
  EXPECT_EQ(3, subsequent_rejections);

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

  outcome = integrator.Solve(problem, adaptive_step_size);
  EXPECT_EQ(TerminationCondition::Done, outcome);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(1e-3 * Metre), Le(2e-3 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(2e-3 * Metre / Second), Le(3e-3 * Metre / Second)));
  EXPECT_EQ(t_initial, solution.back().time.value);
  EXPECT_EQ(steps_backward, solution.size() - steps_forward);
  EXPECT_EQ((1 + initial_rejections) * 4 +
                (steps_backward - 1 + subsequent_rejections) * 3,
            evaluations);
  EXPECT_EQ(1, initial_rejections);
  EXPECT_EQ(11, subsequent_rejections);
}

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest,
       MaxSteps) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      DormandElMikkawyPrince1986RKN434FM<Length>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  Time const period = 2 * π * Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;
  // The number of steps if no step limit is set.
  std::int64_t const steps_forward = 132;

  int evaluations = 0;
  auto const step_size_callback = [](bool tolerable) {};

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
  adaptive_step_size.max_steps = 100;

  auto const outcome = integrator.Solve(problem, adaptive_step_size);
  EXPECT_EQ(TerminationCondition::ReachedMaximalStepCount, outcome);
  EXPECT_THAT(AbsoluteError(
                  x_initial * Cos(ω * (solution.back().time.value - t_initial)),
                      solution.back().positions[0].value),
              AllOf(Ge(8e-4 * Metre), Le(9e-4 * Metre)));
  EXPECT_THAT(AbsoluteError(
                  -v_amplitude *
                      Sin(ω * (solution.back().time.value - t_initial)),
                  solution.back().velocities[0].value),
              AllOf(Ge(1e-3 * Metre / Second), Le(2e-3 * Metre / Second)));
  EXPECT_THAT(solution.back().time.value, Lt(t_final));
  EXPECT_EQ(100, solution.size());

  // Check that a |max_steps| greater than or equal to the unconstrained number
  // of steps has no effect.
  for (std::int64_t const max_steps :
       {steps_forward, steps_forward + 1234}) {
    solution.clear();
    adaptive_step_size.max_steps = steps_forward;
    auto const outcome = integrator.Solve(problem, adaptive_step_size);
    EXPECT_EQ(TerminationCondition::Done, outcome);
    EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
                AllOf(Ge(3e-4 * Metre), Le(4e-4 * Metre)));
    EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
                AllOf(Ge(2e-3 * Metre / Second), Le(3e-3 * Metre / Second)));
    EXPECT_EQ(t_final, solution.back().time.value);
    EXPECT_EQ(steps_forward, solution.size());
  }
}

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest, Singularity) {
  // Integrating the position of an ideal rocket,
  //   x"(t) = m' I_sp / m(t),
  //   x'(0) = 0, x(0) = 0,
  // where m(t) = m₀ - t m'.
  // The solution is
  //   x(t)  = I_sp (t + (t - m₀ / m') log(m₀ / m(t))
  //   x'(t) = I_sp log(m₀ / m(t)) (Циолко́вский's equation).
  // There is a singularity at t = m₀ / m'.
  Variation<Mass> const mass_flow = 1 * Kilogram / Second;
  Mass const initial_mass = 1 * Kilogram;
  SpecificImpulse const specific_impulse = 1 * Newton * Second / Kilogram;
  Instant const t_initial;
  Instant const t_singular = t_initial + initial_mass / mass_flow;
  // After the singularity.
  Instant const t_final = t_initial + 2 * initial_mass / mass_flow;
  auto const mass = [initial_mass, t_initial, mass_flow](Instant const& t) {
    return initial_mass - (t - t_initial) * mass_flow;
  };
  ODE::SystemState initial_state =
      {{0 * Metre}, {0 * Metre / Second}, t_initial};

  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;

  std::vector<ODE::SystemState> solution;
  ODE rocket_equation;
  rocket_equation.compute_acceleration = [&mass, specific_impulse, mass_flow](
      Instant const& t,
      std::vector<Length> const& position,
      not_null<std::vector<Acceleration>*> acceleration) {
    acceleration->back() = mass_flow * specific_impulse / mass(t);
  };
  IntegrationProblem<ODE> problem;
  problem.append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };
  problem.equation = rocket_equation;
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  AdaptiveStepSize<ODE> adaptive_step_size;
  adaptive_step_size.first_time_step = t_final - t_initial;
  adaptive_step_size.safety_factor = 0.9;
  adaptive_step_size.tolerance_to_error_ratio = [length_tolerance,
                                                 speed_tolerance](
      Time const& h, ODE::SystemStateError const& error) {
    return std::min(length_tolerance / Abs(error.position_error[0]),
                    speed_tolerance / Abs(error.velocity_error[0]));
  };

  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      DormandElMikkawyPrince1986RKN434FM<Length>();
  auto const outcome = integrator.Solve(problem, adaptive_step_size);
  EXPECT_EQ(TerminationCondition::VanishingStepSize, outcome);
  EXPECT_EQ(130, solution.size());
  EXPECT_THAT(solution.back().time.value - t_initial,
              AlmostEquals(t_singular - t_initial, 20));
  EXPECT_THAT(solution.back().positions.back().value,
              AlmostEquals(specific_impulse * initial_mass / mass_flow, 711));
}

}  // namespace integrators
}  // namespace principia
