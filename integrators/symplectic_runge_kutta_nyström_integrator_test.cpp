#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Sin;
using quantities::Length;
using quantities::Mass;
using si::Metre;
using si::Milli;
using si::Radian;
using si::Second;
using quantities::Speed;
using quantities::Stiffness;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;
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

class SimpleHarmonicMotionTestInstance {

};

template<typename Integrator>
void Test1000SecondsAt1Millisecond(
    Integrator const& integrator,
    Length const& expected_position_error,
    Speed const& expected_velocity_error) {
  // Long integration, change detector.
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 1000 * Second;
  Time const step = 1 * Milli(Second);
  int const steps = (t_final - t_initial) / step;

  int evaluations = 0;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::SystemState const initial_state = {{q_initial}, {v_initial}, t_initial};
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  problem.append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  integrator.Solve(problem, step);

  EXPECT_EQ(steps, solution.size());
  switch (integrator.composition) {
   case kBA:
   case kABA:
    EXPECT_EQ(steps * integrator.evaluations, evaluations);
    break;
   case kBAB:
    EXPECT_EQ(steps * integrator.evaluations + 1, evaluations);
    break;
  default:
    LOG(FATAL) << "Invalid composition";
  }
  Length q_error;
  Speed v_error;
  for (int i = 1; i <= steps; ++i) {
    Length const q = solution[i - 1].positions[0].value;
    Speed const v = solution[i - 1].velocities[0].value;
    Instant const t = solution[i - 1].time.value;
    EXPECT_THAT(t, AlmostEquals(i * step, 0));
    // TODO(egg): we may need decent trig functions for this sort of thing.
    q_error = std::max(q_error, AbsoluteError(q_initial * Cos(ω * t), q));
    v_error = std::max(v_error, AbsoluteError(v_amplitude * Sin(ω * t), v));
  }
  EXPECT_EQ(expected_position_error, q_error);
  EXPECT_EQ(expected_velocity_error, v_error);
}

}  // namespace

class SymplecticRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(SymplecticRungeKuttaNyströmIntegratorTest,
       HarmonicOscillatorBackAndForth) {
  auto const& aba_integrator = BlanesMoan2002SRKN14A<Length>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;

  int const steps = 100;

  int evaluations = 0;

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

  aba_integrator.Solve(problem, period / 10);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(3E-15 * Metre), Le(4E-15 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(7E-8 * Metre / Second), Le(8E-8 * Metre / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  EXPECT_EQ(steps, solution.size());
  EXPECT_EQ(aba_integrator.evaluations * steps, evaluations);

  evaluations = 0;
  problem.initial_state = &solution.back();
  problem.t_final = t_initial;

  aba_integrator.Solve(problem, -period / 10);
  // The BlanesMoan2002SRKN14A method is reversible: running it backward with
  // the same step size yields the initial state up to roundoff error.
  EXPECT_EQ(x_initial, solution.back().positions[0].value);
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              VanishesBefore(1 * Metre / Second, 2));
  EXPECT_EQ(t_initial, solution.back().time.value);
  EXPECT_EQ(steps, solution.size() - steps);
  EXPECT_EQ(aba_integrator.evaluations * steps, evaluations);

  // Integrate forward with a BAB integrator.
  FixedStepSizeIntegrator<ODE> const& bab_integrator =
      BlanesMoan2002SRKN11B<Length>();
  evaluations = 0;
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  solution.clear();
  bab_integrator.Solve(problem, period / 10);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(4E-14 * Metre), Le(5E-14 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(3E-7 * Metre / Second), Le(4E-7 * Metre / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  EXPECT_EQ(steps, solution.size());
  EXPECT_EQ(aba_integrator.evaluations * steps + 1, evaluations);
}

}  // namespace integrators
}  // namespace principia
