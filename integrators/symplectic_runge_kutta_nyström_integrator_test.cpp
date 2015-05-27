#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using quantities::Acceleration;
using quantities::Length;
using quantities::Mass;
using si::Metre;
using si::Second;
using quantities::Speed;
using quantities::Stiffness;
using testing_utilities::AbsoluteError;
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

}  // namespace

class SymplecticRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(SymplecticRungeKuttaNyströmIntegratorTest,
       HarmonicOscillatorBackAndForth) {
  FixedStepSizeIntegrator<ODE> const& integrator =
      BlanesMoan2002SRKN14A<Length>();
  Length const x_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Time const period = 2 * π * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 10 * period;

  int const steps= 100;

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

  integrator.Solve(problem, period / 10);
  EXPECT_THAT(AbsoluteError(x_initial, solution.back().positions[0].value),
              AllOf(Ge(3E-15 * Metre), Le(4E-15 * Metre)));
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              AllOf(Ge(7E-8 * Metre / Second), Le(8E-8 * Metre / Second)));
  EXPECT_EQ(t_final, solution.back().time.value);
  EXPECT_EQ(steps, solution.size());
  EXPECT_EQ(15 * steps, evaluations);

  evaluations = 0;
  problem.initial_state = &solution.back();
  problem.t_final = t_initial;

  integrator.Solve(problem, -period / 10);
  // The BlanesMoan2002SRKN14A method is reversible: running it backward with
  // the same step size yields the initial state up to roundoff error.
  EXPECT_EQ(x_initial, solution.back().positions[0].value);
  EXPECT_THAT(AbsoluteError(v_initial, solution.back().velocities[0].value),
              VanishesBefore(1 * Metre / Second, 2));
  EXPECT_EQ(t_initial, solution.back().time.value);
  EXPECT_EQ(steps, solution.size() - steps);
  EXPECT_EQ(15 * steps, evaluations);
}

}  // namespace integrators
}  // namespace principia
