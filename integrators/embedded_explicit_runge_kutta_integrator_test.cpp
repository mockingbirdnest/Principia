#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "base/macros.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_integrator {

using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Mass;
using quantities::Sin;
using quantities::SpecificImpulse;
using quantities::Speed;
using quantities::Time;
using quantities::si::Centi;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::ComputeHarmonicOscillatorDerivatives1D;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::StatusIs;
using testing_utilities::operator""_;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::ElementsAreArray;
using ::testing::Lt;
using namespace principia::geometry::_named_quantities;

using ODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length, Speed>;

namespace {

double HarmonicOscillatorToleranceRatio(
    Time const& h,
    ODE::State const& /*state*/,
    ODE::State::Error const& error,
    Length const& q_tolerance,
    Speed const& v_tolerance,
    std::function<void(bool tolerable)> callback) {
  auto const& [position_error, velocity_error] = error;
  double const r = std::min(q_tolerance / Abs(position_error),
                            v_tolerance / Abs(velocity_error));
  callback(r > 1.0);
  return r;
}

}  // namespace

class EmbeddedExplicitRungeKuttaIntegratorTest
    : public ::testing::Test {};

TEST_F(EmbeddedExplicitRungeKuttaIntegratorTest,
       HarmonicOscillatorBackAndForth) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaIntegrator<
          methods::DormandPrince1986RK547FC, ODE>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0.0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;
  int const steps_forward = 50;
  // We integrate backward with double the tolerance.
  int const steps_backward = 46;

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

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_derivative =
      std::bind(ComputeHarmonicOscillatorDerivatives1D,
                _1, _2, _3, &evaluations);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {{x_initial}, {v_initial}}};
  auto const append_state = [&](ODE::State const& state) {
    solution.push_back(state);
  };

  {
    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/t_final - t_initial,
        /*safety_factor=*/0.9);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2, _3,
                  length_tolerance,
                  speed_tolerance,
                  step_size_callback);
    auto instance = integrator.NewInstance(problem,
                                           append_state,
                                           tolerance_to_error_ratio,
                                           parameters);
    auto outcome = instance->Solve(t_final);
    EXPECT_THAT(outcome, StatusIs(termination_condition::Done));
  }
  {
    auto const& [positions, velocities] = solution.back().y;
    EXPECT_THAT(AbsoluteError(x_initial, positions.value),
                IsNear(7.5e-2_(1) * Metre));
    EXPECT_THAT(AbsoluteError(v_initial, velocities.value),
                IsNear(6.8e-2_(1) * Metre / Second));
    EXPECT_EQ(t_final, solution.back().s.value);
    EXPECT_EQ(steps_forward, solution.size());
    EXPECT_EQ(1 + (1 + initial_rejections) * 6 +
                  (steps_forward - 1 + subsequent_rejections) * 6,
              evaluations);
    EXPECT_EQ(1, initial_rejections);
    EXPECT_EQ(1, subsequent_rejections);
  }

  evaluations = 0;
  subsequent_rejections = 0;
  initial_rejections = 0;
  first_step = true;
  problem.initial_state = solution.back();
  {
    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/t_initial - t_final,
        /*safety_factor=*/0.9);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2, _3,
                  2 * length_tolerance,
                  2 * speed_tolerance,
                  step_size_callback);
    auto instance = integrator.NewInstance(
        problem, append_state, tolerance_to_error_ratio, parameters);
    auto outcome = instance->Solve(t_initial);
    EXPECT_THAT(outcome, StatusIs(termination_condition::Done));
  }
  {
    auto const& [positions, velocities] = solution.back().y;
    EXPECT_THAT(AbsoluteError(x_initial, positions.value),
                IsNear(2.1e-1_(1) * Metre));
    EXPECT_THAT(AbsoluteError(v_initial, velocities.value),
                IsNear(6.8e-2_(1) * Metre / Second));
    EXPECT_EQ(t_initial, solution.back().s.value);
    EXPECT_EQ(steps_backward, solution.size() - steps_forward);
    EXPECT_EQ(1 + (1 + initial_rejections) * 6 +
                  (steps_backward - 1 + subsequent_rejections) * 6,
              evaluations);
    EXPECT_EQ(1, initial_rejections);
    EXPECT_EQ(1, subsequent_rejections);
  }
}

TEST_F(EmbeddedExplicitRungeKuttaIntegratorTest, MaxSteps) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaIntegrator<
          methods::DormandPrince1986RK547FC, ODE>();
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
  std::int64_t const steps_forward = 50;

  auto const step_size_callback = [](bool tolerable) {};

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_derivative =
      std::bind(ComputeHarmonicOscillatorDerivatives1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {{x_initial}, {v_initial}}};
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };
  AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9,
      /*max_steps=*/40,
      /*last_step_is_exact=*/true);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2, _3,
                length_tolerance,
                speed_tolerance,
                step_size_callback);

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  auto const outcome = instance->Solve(t_final);

  auto const& [positions, velocities] = solution.back().y;
  EXPECT_THAT(outcome,
              StatusIs(termination_condition::ReachedMaximalStepCount));
  EXPECT_THAT(AbsoluteError(
                  x_initial * Cos(ω * (solution.back().s.value - t_initial)),
                  positions.value),
              IsNear(5.9e-2_(1) * Metre));
  EXPECT_THAT(AbsoluteError(
                  -v_amplitude *
                      Sin(ω * (solution.back().s.value - t_initial)),
                  velocities.value),
              IsNear(5.6e-2_(1) * Metre / Second));
  EXPECT_THAT(solution.back().s.value, Lt(t_final));
  EXPECT_EQ(40, solution.size());

  // Check that a |max_steps| greater than or equal to the unconstrained number
  // of steps has no effect.
  for (std::int64_t const max_steps :
       {steps_forward, steps_forward + 1234}) {
    solution.clear();
    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/t_final - t_initial,
        /*safety_factor=*/0.9,
        /*max_steps=*/max_steps,
        /*last_step_is_exact=*/true);
    auto const instance = integrator.NewInstance(problem,
                                                 append_state,
                                                 tolerance_to_error_ratio,
                                                 parameters);
    auto const outcome = instance->Solve(t_final);
    auto const& [positions, velocities] = solution.back().y;
    EXPECT_THAT(outcome, StatusIs(termination_condition::Done));
    EXPECT_THAT(AbsoluteError(x_initial, positions.value),
                IsNear(7.5e-2_(1) * Metre));
    EXPECT_THAT(AbsoluteError(v_initial, velocities.value),
                IsNear(6.8e-2_(1) * Metre / Second));
    EXPECT_EQ(t_final, solution.back().s.value);
    EXPECT_EQ(steps_forward, solution.size());
  }
}

TEST_F(EmbeddedExplicitRungeKuttaIntegratorTest, Singularity) {
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

  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;

  std::vector<ODE::State> solution;
  ODE rocket_equation;
  rocket_equation.compute_derivative = [&mass, specific_impulse, mass_flow](
      Instant const& t,
      ODE::DependentVariables const& dependent_variables,
      ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
    auto const& [q, v] = dependent_variables;
    auto& [qʹ, vʹ] = dependent_variable_derivatives;
    qʹ = v;
    vʹ = mass_flow * specific_impulse / mass(t);
    return absl::OkStatus();
  };
  InitialValueProblem<ODE> problem;
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };
  problem.equation = rocket_equation;
  problem.initial_state = {t_initial, {{0 * Metre}, {0 * Metre / Second}}};
  AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio = [length_tolerance, speed_tolerance](
      Time const& h,
      ODE::State const& /*state*/,
      ODE::State::Error const& error) {
    auto const& [position_error, velocity_error] = error;
    return std::min(length_tolerance / Abs(position_error),
                    speed_tolerance / Abs(velocity_error));
  };

  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaIntegrator<
          methods::DormandPrince1986RK547FC, ODE>();

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  auto const outcome = instance->Solve(t_final);

  auto const& [position, velocity] = solution.back().y;
  EXPECT_THAT(outcome, StatusIs(termination_condition::VanishingStepSize));
  EXPECT_EQ(101, solution.size());
  EXPECT_THAT(solution.back().s.value - t_initial,
              AlmostEquals(t_singular - t_initial, 4));
  EXPECT_THAT(position.value,
              AlmostEquals(specific_impulse * initial_mass / mass_flow, 155));
}

TEST_F(EmbeddedExplicitRungeKuttaIntegratorTest, Restart) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaIntegrator<
          methods::DormandPrince1986RK547FC, ODE>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Time const duration = 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;

  auto const step_size_callback = [](bool tolerable) {};

  std::vector<ODE::State> solution1;
  {
    ODE harmonic_oscillator;
    harmonic_oscillator.compute_derivative =
        std::bind(ComputeHarmonicOscillatorDerivatives1D,
                  _1, _2, _3, /*evaluations=*/nullptr);
    InitialValueProblem<ODE> problem;
    problem.equation = harmonic_oscillator;
    problem.initial_state = {t_initial, {{x_initial}, {v_initial}}};
    auto const append_state = [&solution1](ODE::State const& state) {
      solution1.push_back(state);
    };

    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/duration,
        /*safety_factor=*/0.9,
        /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
        /*last_step_is_exact=*/false);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2, _3,
                  length_tolerance,
                  speed_tolerance,
                  step_size_callback);

    auto const instance = integrator.NewInstance(problem,
                                                 append_state,
                                                 tolerance_to_error_ratio,
                                                 parameters);
    auto outcome = instance->Solve(t_initial + duration);
    EXPECT_THAT(outcome, StatusIs(termination_condition::Done));

    // Check that the time step has been updated.
    EXPECT_EQ(49, solution1.size());
    EXPECT_THAT(
        solution1[solution1.size() - 1].s.value -
            solution1[solution1.size() - 2].s.value,
        AlmostEquals(1.247'142'266'261'910'49 * Second, 0));

    // Restart the integration.
    outcome = instance->Solve(t_initial + 2.0 * duration);
    EXPECT_THAT(outcome, StatusIs(termination_condition::Done));

    // Check that the time step has been updated again.
    EXPECT_EQ(99, solution1.size());
    EXPECT_THAT(
        solution1[solution1.size() - 1].s.value -
            solution1[solution1.size() - 2].s.value,
        AlmostEquals(1.237'882'807'009'299'31 * Second, 0));
  }

  // Do it again in one call to |Solve| and check associativity.
  std::vector<ODE::State> solution2;
  {
    ODE harmonic_oscillator;
    harmonic_oscillator.compute_derivative =
        std::bind(ComputeHarmonicOscillatorDerivatives1D,
                  _1, _2, _3, /*evaluations=*/nullptr);
    InitialValueProblem<ODE> problem;
    problem.equation = harmonic_oscillator;
    problem.initial_state = {t_initial, {{x_initial}, {v_initial}}};
    auto const append_state = [&solution2](ODE::State const& state) {
      solution2.push_back(state);
    };

    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/duration,
        /*safety_factor=*/0.9,
        /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
        /*last_step_is_exact=*/false);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2, _3,
                  length_tolerance,
                  speed_tolerance,
                  step_size_callback);

    auto const instance = integrator.NewInstance(problem,
                                                 append_state,
                                                 tolerance_to_error_ratio,
                                                 parameters);
    auto outcome = instance->Solve(t_initial + 2.0 * duration);
    EXPECT_THAT(outcome, StatusIs(termination_condition::Done));
  }

  EXPECT_THAT(solution2, ElementsAreArray(solution1));
}

#if 0
TEST_F(EmbeddedExplicitRungeKuttaIntegratorTest, Serialization) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaIntegrator<
          methods::DormandPrince1986RK547FC, Length, Speed>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;

  auto const step_size_callback = [](bool tolerable) {};

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_derivative =
      std::bind(ComputeHarmonicOscillatorDerivatives1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{{x_initial}, {v_initial}}, t_initial};
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };
  AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2,
                length_tolerance,
                speed_tolerance,
                step_size_callback);

  auto const instance1 = integrator.NewInstance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters);
  serialization::IntegratorInstance message1;
  instance1->WriteToMessage(&message1);
  auto const instance2 =
      AdaptiveStepSizeIntegrator<ODE>::Instance::ReadFromMessage(
          message1,
          harmonic_oscillator,
          append_state,
          tolerance_to_error_ratio);
  serialization::IntegratorInstance message2;
  instance2->WriteToMessage(&message2);
  EXPECT_THAT(message1, EqualsProto(message2));
}
#endif

}  // namespace internal_embedded_explicit_runge_kutta_integrator

// Reopen this namespace to allow printing out the system state.
namespace internal_ordinary_differential_equations {

void PrintTo(
    typename internal_embedded_explicit_runge_kutta_integrator::ODE::
        State const& state,
    std::ostream* const out) {
  auto const& [position, velocity] = state.y;
  *out << "\nTime: " << state.s << "\n";
  *out << "Position: " << position << "\n";
  *out << "Velocity: " << velocity << "\n";
}

}  // namespace internal_ordinary_differential_equations

}  // namespace integrators
}  // namespace principia
