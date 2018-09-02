
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_nyström_integrator {

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
using testing_utilities::ComputeHarmonicOscillatorAcceleration1D;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::ElementsAreArray;
using ::testing::Lt;

using ODE = SpecialSecondOrderDifferentialEquation<Length>;

namespace {

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
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::DormandالمكاوىPrince1986RKN434FM,
          Length>();
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
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{x_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  {
    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/t_final - t_initial,
        /*safety_factor=*/0.9);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2,
                  length_tolerance,
                  speed_tolerance,
                  step_size_callback);
    auto instance = integrator.NewInstance(problem,
                                           append_state,
                                           tolerance_to_error_ratio,
                                           parameters);
    auto outcome = instance->Solve(t_final);
    EXPECT_EQ(termination_condition::Done, outcome.error());
  }
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              IsNear(3.5e-4 * Metre));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              IsNear(2.8e-3 * Metre / Second));
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
  problem.initial_state = solution.back();
  {
    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/t_initial - t_final,
        /*safety_factor=*/0.9);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2,
                  2 * length_tolerance,
                  2 * speed_tolerance,
                  step_size_callback);
    auto instance = integrator.NewInstance(
        problem, append_state, tolerance_to_error_ratio, parameters);
    auto outcome = instance->Solve(t_initial);
    EXPECT_EQ(termination_condition::Done, outcome.error());
  }
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              IsNear(1.2e-3 * Metre));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              IsNear(2.5e-3 * Metre / Second));
  EXPECT_EQ(t_initial, solution.back().time.value);
  EXPECT_EQ(steps_backward, solution.size() - steps_forward);
  EXPECT_EQ((1 + initial_rejections) * 4 +
                (steps_backward - 1 + subsequent_rejections) * 3,
            evaluations);
  EXPECT_EQ(1, initial_rejections);
  EXPECT_EQ(11, subsequent_rejections);
}

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest, MaxSteps) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::DormandالمكاوىPrince1986RKN434FM,
          Length>();
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

  auto const step_size_callback = [](bool tolerable) {};

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{x_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };
  AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9,
      /*max_steps=*/100,
      /*last_step_is_exact=*/true);
  auto const tolerance_to_error_ratio =
      std::bind(HarmonicOscillatorToleranceRatio,
                _1, _2,
                length_tolerance,
                speed_tolerance,
                step_size_callback);

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  auto const outcome = instance->Solve(t_final);

  EXPECT_EQ(termination_condition::ReachedMaximalStepCount, outcome.error());
  EXPECT_THAT(AbsoluteError(
                  x_initial * Cos(ω * (solution.back().time.value - t_initial)),
                      solution.back().positions[0].value),
              IsNear(9.0e-4 * Metre));
  EXPECT_THAT(AbsoluteError(
                  -v_amplitude *
                      Sin(ω * (solution.back().time.value - t_initial)),
                  solution.back().velocities[0].value),
              IsNear(1.9e-3 * Metre / Second));
  EXPECT_THAT(solution.back().time.value, Lt(t_final));
  EXPECT_EQ(100, solution.size());

  // Check that a |max_steps| greater than or equal to the unconstrained number
  // of steps has no effect.
  for (std::int64_t const max_steps :
       {steps_forward, steps_forward + 1234}) {
    solution.clear();
    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/t_final - t_initial,
        /*safety_factor=*/0.9,
        /*max_steps=*/steps_forward,
        /*last_step_is_exact=*/true);
    auto const instance = integrator.NewInstance(problem,
                                                 append_state,
                                                 tolerance_to_error_ratio,
                                                 parameters);
    auto const outcome = instance->Solve(t_final);
    EXPECT_EQ(termination_condition::Done, outcome.error());
    EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
                IsNear(3.6e-4 * Metre));
    EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
                IsNear(2.8e-3 * Metre / Second));
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

  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;

  std::vector<ODE::SystemState> solution;
  ODE rocket_equation;
  rocket_equation.compute_acceleration = [&mass, specific_impulse, mass_flow](
      Instant const& t,
      std::vector<Length> const& position,
      std::vector<Acceleration>& acceleration) {
    acceleration.back() = mass_flow * specific_impulse / mass(t);
    return Status::OK;
  };
  IntegrationProblem<ODE> problem;
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };
  problem.equation = rocket_equation;
  problem.initial_state = {{0 * Metre}, {0 * Metre / Second}, t_initial};
  AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio = [length_tolerance, speed_tolerance](
      Time const& h, ODE::SystemStateError const& error) {
    return std::min(length_tolerance / Abs(error.position_error[0]),
                    speed_tolerance / Abs(error.velocity_error[0]));
  };

  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::DormandالمكاوىPrince1986RKN434FM,
          Length>();

  auto const instance = integrator.NewInstance(problem,
                                               append_state,
                                               tolerance_to_error_ratio,
                                               parameters);
  auto const outcome = instance->Solve(t_final);

  EXPECT_EQ(termination_condition::VanishingStepSize, outcome.error());
  EXPECT_EQ(132, solution.size());
  EXPECT_THAT(solution.back().time.value - t_initial,
              AlmostEquals(t_singular - t_initial, 15));
  EXPECT_THAT(solution.back().positions.back().value,
              AlmostEquals(specific_impulse * initial_mass / mass_flow, 540));
}

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest, Restart) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::DormandالمكاوىPrince1986RKN434FM,
          Length>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  Time const period = 2 * π * Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
  Time const duration = 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;
  // The number of steps if no step limit is set.
  std::int64_t const steps_forward = 132;

  auto const step_size_callback = [](bool tolerable) {};

  std::vector<ODE::SystemState> solution1;
  {
    ODE harmonic_oscillator;
    harmonic_oscillator.compute_acceleration =
        std::bind(ComputeHarmonicOscillatorAcceleration1D,
                  _1, _2, _3, /*evaluations=*/nullptr);
    IntegrationProblem<ODE> problem;
    problem.equation = harmonic_oscillator;
    problem.initial_state = {{x_initial}, {v_initial}, t_initial};
    auto const append_state = [&solution1](ODE::SystemState const& state) {
      solution1.push_back(state);
    };

    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/duration,
        /*safety_factor=*/0.9,
        /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
        /*last_step_is_exact=*/false);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2,
                  length_tolerance,
                  speed_tolerance,
                  step_size_callback);

    auto const instance = integrator.NewInstance(problem,
                                                 append_state,
                                                 tolerance_to_error_ratio,
                                                 parameters);
    auto outcome = instance->Solve(t_initial + duration);
    EXPECT_EQ(termination_condition::Done, outcome.error());

    // Check that the time step has been updated.
    EXPECT_EQ(131, solution1.size());
    EXPECT_THAT(
        solution1[solution1.size() - 1].time.value -
            solution1[solution1.size() - 2].time.value,
        AlmostEquals(0.509'363'975'335'290'320 * Second, 0));

    // Restart the integration.
    outcome = instance->Solve(t_initial + 2.0 * duration);
    EXPECT_EQ(termination_condition::Done, outcome.error());

    // Check that the time step has been updated again.
    EXPECT_EQ(261, solution1.size());
    EXPECT_THAT(
        solution1[solution1.size() - 1].time.value -
            solution1[solution1.size() - 2].time.value,
        AlmostEquals(0.506'410'259'195'249'068 * Second, 0));
  }

  // Do it again in one call to |Solve| and check associativity.
  std::vector<ODE::SystemState> solution2;
  {
    ODE harmonic_oscillator;
    harmonic_oscillator.compute_acceleration =
        std::bind(ComputeHarmonicOscillatorAcceleration1D,
                  _1, _2, _3, /*evaluations=*/nullptr);
    IntegrationProblem<ODE> problem;
    problem.equation = harmonic_oscillator;
    problem.initial_state = {{x_initial}, {v_initial}, t_initial};
    auto const append_state = [&solution2](ODE::SystemState const& state) {
      solution2.push_back(state);
    };

    AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
        /*first_time_step=*/duration,
        /*safety_factor=*/0.9,
        /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
        /*last_step_is_exact=*/false);
    auto const tolerance_to_error_ratio =
        std::bind(HarmonicOscillatorToleranceRatio,
                  _1, _2,
                  length_tolerance,
                  speed_tolerance,
                  step_size_callback);

    auto const instance = integrator.NewInstance(problem,
                                                 append_state,
                                                 tolerance_to_error_ratio,
                                                 parameters);
    auto outcome = instance->Solve(t_initial + 2.0 * duration);
    EXPECT_EQ(termination_condition::Done, outcome.error());
  }

  EXPECT_THAT(solution2, ElementsAreArray(solution1));
}

TEST_F(EmbeddedExplicitRungeKuttaNyströmIntegratorTest, Serialization) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          methods::DormandالمكاوىPrince1986RKN434FM,
          Length>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;
  Length const length_tolerance = 1 * Milli(Metre);
  Speed const speed_tolerance = 1 * Milli(Metre) / Second;

  auto const step_size_callback = [](bool tolerable) {};

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{x_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
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

}  // namespace internal_embedded_explicit_runge_kutta_nyström_integrator

// Reopen this namespace to allow printing out the system state.
namespace internal_ordinary_differential_equations {

void PrintTo(
    typename internal_embedded_explicit_runge_kutta_nyström_integrator::ODE::
        SystemState const& system_state,
    std::ostream* const out) {
  *out << "\nTime: " << system_state.time << "\n";
  *out << "Positions:\n";
  for (int i = 0; i < system_state.positions.size(); ++i) {
    *out << "  " << i << ": " << system_state.positions[i] << "\n";
  }
  *out << "Velocities:\n";
  for (int i = 0; i < system_state.velocities.size(); ++i) {
    *out << "  " << i << ": " << system_state.velocities[i] << "\n";
  }
}

}  // namespace internal_ordinary_differential_equations

}  // namespace integrators
}  // namespace principia
