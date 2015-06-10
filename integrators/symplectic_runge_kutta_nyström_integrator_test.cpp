#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <vector>

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
using quantities::Sin;
using quantities::Length;
using quantities::Mass;
using quantities::Pow;
using quantities::Power;
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

#define INSTANCE(integrator,                                    \
                 beginning_of_convergence,                      \
                 expected_position_error,                       \
                 expected_velocity_error,                       \
                 expected_energy_error)                         \
    SimpleHarmonicMotionTestInstance(integrator<Length>(),      \
                                     #integrator,               \
                                     beginning_of_convergence,  \
                                     expected_position_error,   \
                                     expected_velocity_error,   \
                                     expected_energy_error)

namespace integrators {

using ODE = SpecialSecondOrderDifferentialEquation<Length>;

namespace {

// TODO(egg): use the one from testing_utilities/integration again when everyone
// uses |Instant|s.
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
      name_(name) {};

  std::string const& name() const {
    return name_;
  };

  void Run1000SecondsAt1Millisecond() const {
    test_1000_seconds_at_1_millisecond_();
  };

  void RunConvergence() const {
    test_convergence_();
  };

  void RunSymplecticity()  const {
    test_symplecticity_();
  };

 private:
  std::function<void()> test_1000_seconds_at_1_millisecond_;
  std::function<void()> test_convergence_;
  std::function<void()> test_symplecticity_;
  std::string name_;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         SimpleHarmonicMotionTestInstance instance) {
  return stream << instance.name();
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
  Time const period = 2 * π * Second;
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
  for (int i = 1; i <= steps; ++i) {
    Length const q = solution[i - 1].positions[0].value;
    Speed const v = solution[i - 1].velocities[0].value;
    Time const t = solution[i - 1].time.value - t_initial;
    EXPECT_THAT(t, AlmostEquals(i * step, 0));
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
  Time const period = 2 * π * Second;
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
                _1, _2, _3, nullptr);
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
     RelativeError(((integrator.order + 1) / 2) * 2,
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
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Time const period = 2 * π * Second;
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
                _1, _2, _3, nullptr);
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

std::vector<SimpleHarmonicMotionTestInstance> Instances() {
  return {INSTANCE(McLachlanAtela1992Order4Optimal,
                   1.0 * Second,
                   +1.88161985992252310e-13 * Metre,
                   +1.88491583452687910e-13 * Metre / Second,
                   +7.52285331973023830e-07 * Joule),
          INSTANCE(McLachlanAtela1992Order5Optimal,
                   1.1 * Second,
                   +7.51005160837259210e-14 * Metre,
                   +7.50823014872281650e-14 * Metre / Second,
                   +3.06327349042234690e-08 * Joule)};
}

}  // namespace

class SymplecticRungeKuttaNyströmIntegratorTest
    : public ::testing::TestWithParam<SimpleHarmonicMotionTestInstance> {};

INSTANTIATE_TEST_CASE_P(SymplecticRungeKuttaNyströmIntegratorTests,
                        SymplecticRungeKuttaNyströmIntegratorTest,
                        ValuesIn(Instances()));

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Symplecticity) {
  LOG(INFO) << GetParam();
  GetParam().RunSymplecticity();
}

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Convergence) {
  LOG(INFO) << GetParam();
  GetParam().RunConvergence();
}

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, LongIntegration) {
  LOG(INFO) << GetParam();
  GetParam().Run1000SecondsAt1Millisecond();
}

}  // namespace integrators
}  // namespace principia
