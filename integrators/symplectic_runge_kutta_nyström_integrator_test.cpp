﻿#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <vector>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Energy;
using quantities::Length;
using quantities::Mass;
using quantities::Pow;
using quantities::Power;
using quantities::Sin;
using si::Joule;
using si::Kilogram;
using si::Metre;
using si::Milli;
using si::Radian;
using si::Second;
using quantities::Speed;
using quantities::Stiffness;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;
using testing_utilities::VanishesBefore;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;
using ::testing::ValuesIn;

#define INSTANCE(integrator,                                      \
                 beginning_of_convergence,                        \
                 expected_position_error,                         \
                 expected_velocity_error,                         \
                 expected_energy_error)                           \
    SimpleHarmonicMotionTestInstance(integrator<Length>(),        \
                                     #integrator,                 \
                                     (beginning_of_convergence),  \
                                     (expected_position_error),   \
                                     (expected_velocity_error),   \
                                     (expected_energy_error))

namespace integrators {

using ODE = SpecialSecondOrderDifferentialEquation<Length>;

namespace {

// TODO(egg): use the one from testing_utilities/integration again when everyone
// uses |Instant|s.
// Increments |*evaluations| if |evaluations| is not null.
void ComputeHarmonicOscillatorAcceleration(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* const result,
    int* evaluations) {
  (*result)[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
  if (evaluations != nullptr) {
    ++*evaluations;
  }
}

// Long integration, change detector.  Also tests the number of steps, their
// spacing, and the number of evaluations.
template<typename Integrator>
void Test1000SecondsAt1Millisecond(
    Integrator const& integrator,
    Length const& expected_position_error,
    Speed const& expected_velocity_error) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 1000 * Second;
  Time const step = 1 * Milli(Second);
  int const steps = static_cast<int>((t_final - t_initial) / step) - 1;

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
  for (int i = 0; i < steps; ++i) {
    Length const q = solution[i].positions[0].value;
    Speed const v = solution[i].velocities[0].value;
    Time const t = solution[i].time.value - t_initial;
    EXPECT_THAT(t, AlmostEquals((i + 1) * step, 0));
    // TODO(egg): we may need decent trig functions for this sort of thing.
    q_error = std::max(q_error, AbsoluteError(q_initial * Cos(ω * t), q));
    v_error = std::max(v_error, AbsoluteError(-v_amplitude * Sin(ω * t), v));
  }
  EXPECT_EQ(expected_position_error, q_error);
  EXPECT_EQ(expected_velocity_error, v_error);
}

// Integrates with diminishing step sizes, and checks the order of convergence.
template<typename Integrator>
void TestConvergence(Integrator const& integrator,
                     Time const& beginning_of_convergence) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 100 * Second;

  Time step = beginning_of_convergence;
  int const step_sizes = 50;
  double const step_reduction = 1.1;
  std::vector<double> log_step_sizes;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_q_errors;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_p_errors;
  log_step_sizes.reserve(step_sizes);

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, nullptr /*evaluations*/);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::SystemState const initial_state = {{q_initial}, {v_initial}, t_initial};
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  ODE::SystemState final_state;
  problem.append_state = [&final_state](ODE::SystemState const& state) {
    final_state = state;
  };

  for (int i = 0; i < step_sizes; ++i, step /= step_reduction) {
    integrator.Solve(problem, step);
    Time const t = final_state.time.value - t_initial;
    Length const& q = final_state.positions[0].value;
    Speed const& v = final_state.velocities[0].value;
    double const log_q_error = std::log10(
        AbsoluteError(q / q_initial, Cos(ω * t)));
    double const log_p_error = std::log10(
        AbsoluteError(v / v_amplitude, -Sin(ω * t)));
    if (log_q_error <= -13 || log_p_error <= -13) {
      // If we keep going the effects of finite precision will drown out
      // convergence.
      break;
    }
    log_step_sizes.push_back(std::log10(step / Second));
    log_q_errors.push_back(log_q_error);
    log_p_errors.push_back(log_p_error);
  }
  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;

  EXPECT_THAT(RelativeError(integrator.order, q_convergence_order),
              Lt(0.02));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.99), Lt(1.01)));
  double const v_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_order;
  LOG(INFO) << "Correlation            : " << v_correlation;
  // SPRKs with odd convergence order have a higher convergence order in p.
  EXPECT_THAT(
      RelativeError(integrator.order + (integrator.order % 2),
                    v_convergence_order),
      Lt(0.02));
  EXPECT_THAT(v_correlation, AllOf(Gt(0.99), Lt(1.01)));
}

// Test that the error in energy does not correlate with the number of steps
// taken, and checks the actual value of the error.
template<typename Integrator>
void TestSymplecticity(Integrator const& integrator,
                       Energy const& expected_energy_error) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 500 * Second;
  Time const step = 0.2 * Second;

  Mass const m = 1 * Kilogram;
  Stiffness const k = SIUnit<Stiffness>();
  Energy const initial_energy =
      0.5 * m * Pow<2>(v_initial) + 0.5 * k * Pow<2>(q_initial);

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, nullptr /*evaluations*/);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::SystemState const initial_state = {{q_initial}, {v_initial}, t_initial};
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  problem.append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  integrator.Solve(problem, step);

  std::size_t const length = solution.size();
  std::vector<Energy> energy_error(length);
  std::vector<Time> time(length);
  Energy max_energy_error;
  for (std::size_t i = 0; i < length; ++i) {
    Length const q_i   = solution[i].positions[0].value;
    Speed const v_i = solution[i].velocities[0].value;
    time[i] = solution[i].time.value - t_initial;
    energy_error[i] =
        AbsoluteError(initial_energy,
                      0.5 * m * Pow<2>(v_i) + 0.5 * k * Pow<2>(q_i));
    max_energy_error = std::max(energy_error[i], max_energy_error);
  }
  double const correlation =
      PearsonProductMomentCorrelationCoefficient(time, energy_error);
  LOG(INFO) << "Correlation between time and energy error : " << correlation;
  EXPECT_THAT(correlation, Lt(2E-3));
  Power const slope = Slope(time, energy_error);
  LOG(INFO) << "Slope                                     : " << slope;
  EXPECT_THAT(Abs(slope), Lt(2E-6 * SIUnit<Power>()));
  LOG(INFO) << "Maximum energy error                      : " <<
      max_energy_error;
  EXPECT_EQ(expected_energy_error, max_energy_error);
}

// If the integrator is time-reversible, checks that integrating back and
// forth with a large time step yields the initial value.  If it is not, checks
// that there is an error when integrating back and forth.
template<typename Integrator>
void TestTimeReversibility(Integrator const& integrator) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 100 * Second;
  Time const step = 1 * Second;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
      std::bind(ComputeHarmonicOscillatorAcceleration,
                _1, _2, _3, nullptr /*evaluations*/);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::SystemState const initial_state = {{q_initial}, {v_initial}, t_initial};
  ODE::SystemState final_state;
  problem.initial_state = &initial_state;
  problem.t_final = t_final;
  problem.append_state = [&final_state](ODE::SystemState const& state) {
    final_state = state;
  };

  integrator.Solve(problem, step);

  problem.initial_state = &final_state;
  problem.t_final = t_initial;

  integrator.Solve(problem, -step);

  EXPECT_EQ(t_initial, final_state.time.value);
  if (integrator.time_reversible) {
    EXPECT_THAT(final_state.positions[0].value,
                AlmostEquals(q_initial, 0, 8));
    EXPECT_THAT(final_state.velocities[0].value,
                VanishesBefore(v_amplitude, 0, 16));
  } else {
    EXPECT_THAT(AbsoluteError(q_initial,
                              final_state.positions[0].value),
                Gt(1E-4 * Metre));
    EXPECT_THAT(AbsoluteError(v_initial,
                              final_state.velocities[0].value),
                Gt(1E-4 * Metre / Second));
  }
}

struct SimpleHarmonicMotionTestInstance {
 public:
  template<typename Integrator>
  SimpleHarmonicMotionTestInstance(Integrator const& integrator,
                                   std::string const& name,
                                   Time const& beginning_of_convergence,
                                   Length const& expected_position_error,
                                   Speed const& expected_velocity_error,
                                   Energy const& expected_energy_error)
      : test_1000_seconds_at_1_millisecond_(
            std::bind(Test1000SecondsAt1Millisecond<Integrator>,
                      integrator,
                      expected_position_error,
                      expected_velocity_error)),
        test_convergence_(
            std::bind(TestConvergence<Integrator>,
                      integrator,
                      beginning_of_convergence)),
        test_symplecticity_(
            std::bind(TestSymplecticity<Integrator>,
                      integrator,
                      expected_energy_error)),
        test_time_reversibility_(
            std::bind(TestTimeReversibility<Integrator>,
                      integrator)),
        name_(name) {}

  std::string const& name() const {
    return name_;
  }

  void Run1000SecondsAt1Millisecond() const {
    test_1000_seconds_at_1_millisecond_();
  }

  void RunConvergence() const {
    test_convergence_();
  }

  void RunSymplecticity()  const {
    test_symplecticity_();
  }

  void RunTimeReversibility() const {
    test_time_reversibility_();
  }

 private:
  std::function<void()> test_1000_seconds_at_1_millisecond_;
  std::function<void()> test_convergence_;
  std::function<void()> test_symplecticity_;
  std::function<void()> test_time_reversibility_;
  std::string name_;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         SimpleHarmonicMotionTestInstance instance) {
  return stream << instance.name();
}

std::vector<SimpleHarmonicMotionTestInstance> Instances() {
  return {INSTANCE(McLachlanAtela1992Order4Optimal,
                   1.0 * Second,
                   +1.88161985992252310e-13 * Metre,
                   +1.88491583452687910e-13 * Metre / Second,
                   +7.52285331973023830e-07 * Joule),
          INSTANCE(McLachlan1995SB3A4,
                   1.0 * Second,
                   +1.21597176772070270e-13 * Metre,
                   +1.21782792183999790e-13 * Metre / Second,
                   +2.52639347009253610e-06 * Joule),
          INSTANCE(McLachlan1995SB3A5,
                   1.0 * Second,
                   +1.37754391227318250e-13 * Metre,
                   +1.37848066295021000e-13 * Metre / Second,
                   +1.70551544109720510e-07 * Joule),
          INSTANCE(BlanesMoan2002SRKN6B,
                   1.0 * Second,
                   +1.18405285576272950e-13 * Metre,
                   +1.18564880136062810e-13 * Metre / Second,
                   +1.55706381121945010e-09 * Joule),
          INSTANCE(McLachlanAtela1992Order5Optimal,
                   1.1 * Second,
                   +7.51005160837259210e-14 * Metre,
                   +7.50823014872281650e-14 * Metre / Second,
                   +3.06327349042234690e-08 * Joule),
          INSTANCE(OkunborSkeel1994Order6Method13,
                   1.1 * Second,
                   +1.54723456269323380e-13 * Metre,
                   +1.54959378662056220e-13 * Metre / Second,
                   +2.28626773068896230e-09 * Joule),
          INSTANCE(BlanesMoan2002SRKN11B,
                   1.0 * Second,
                   +1.11972930927350940e-13 * Metre,
                   +1.12035380972486110e-13 * Metre / Second,
                   +9.14945896823837760e-12 * Joule),
          INSTANCE(BlanesMoan2002SRKN14A,
                   1.0 * Second,
                   +1.11001485780803930e-13 * Metre,
                   +1.11063935825939100e-13 * Metre / Second,
                   +6.29052365752613700e-13 * Joule)};
}

}  // namespace

// Beware!  Unicode minefield ahead!
// The name of this class and subsequent tests is not what you think.  It turns
// out that if we use a 'LATIN SMALL LETTER O WITH DIAERESIS' (U+00F6) some of
// the tools that process tests interpret it as a 'DIVISION SIGN' (U+00F7) and
// fail.  To avoid these problems, we use a 'LATIN SMALL LETTER O' (U+006F)
// followed by a 'COMBINING DIAERESIS' (U+0308) as it seems that the code points
// in the "upper half" are those that cause trouble.
class SymplecticRungeKuttaNyströmIntegratorTest
    : public ::testing::TestWithParam<SimpleHarmonicMotionTestInstance> {
 public:
  SymplecticRungeKuttaNyströmIntegratorTest() {
    google::LogToStderr();
  }
};

INSTANTIATE_TEST_CASE_P(SymplecticRungeKuttaNyströmIntegratorTests,
                        SymplecticRungeKuttaNyströmIntegratorTest,
                        ValuesIn(Instances()));

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, TimeReversibility) {
  LOG(INFO) << GetParam();
  GetParam().RunTimeReversibility();
}

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Symplecticity) {
  LOG(INFO) << GetParam();
  GetParam().RunSymplecticity();
}

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Convergence) {
  LOG(INFO) << GetParam();
  GetParam().RunConvergence();
}

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, LongIntegration) {
  LOG(INFO) << GetParam();
  GetParam().Run1000SecondsAt1Millisecond();
}

}  // namespace integrators
}  // namespace principia
