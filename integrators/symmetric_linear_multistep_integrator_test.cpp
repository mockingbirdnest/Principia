#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include <algorithm>
#include <filesystem>
#include <vector>
#include <string>

#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/integrators.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "mathematica/logger.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace integrators {

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;
using ::testing::ValuesIn;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::mathematica::_logger;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_integration;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_statistics;
using namespace principia::testing_utilities::_vanishes_before;

#define INSTANCE(integrator,                                          \
                 beginning_of_convergence,                            \
                 expected_position_error,                             \
                 expected_velocity_error,                             \
                 expected_energy_error)                               \
  SimpleHarmonicMotionTestInstance(                                   \
      SymmetricLinearMultistepIntegrator<methods::integrator, ODE>(), \
      #integrator,                                                    \
      (beginning_of_convergence),                                     \
      (expected_position_error),                                      \
      (expected_velocity_error),                                      \
      (expected_energy_error))

using ODE = SpecialSecondOrderDifferentialEquation<Length>;

namespace {

struct SimpleHarmonicMotionTestInstance final {
  template<typename Integrator>
  SimpleHarmonicMotionTestInstance(Integrator const& integrator,
                                   std::string const& name,
                                   Time const& beginning_of_convergence,
                                   Length const& expected_position_error,
                                   Speed const& expected_velocity_error,
                                   Energy const& expected_energy_error)
      : integrator(integrator),
        order(integrator.order),
        name(name),
        beginning_of_convergence(beginning_of_convergence),
        expected_position_error(expected_position_error),
        expected_velocity_error(expected_velocity_error),
        expected_energy_error(expected_energy_error) {}

  FixedStepSizeIntegrator<ODE> const& integrator;
  int const order;
  std::string const name;
  Time const beginning_of_convergence;
  Length const expected_position_error;
  Speed const expected_velocity_error;
  Energy const expected_energy_error;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         SimpleHarmonicMotionTestInstance const& instance) {
  return stream << instance.name;
}

// Not testing QuinlanTremaine1990Order14 as its characteristics cannot even be
// computed.
std::vector<SimpleHarmonicMotionTestInstance> Instances() {
  // The `beginning_of_convergence` below were carefully chosen using
  // Mathematica to only select the domain where `p` and `step` are properly
  // correlated.
  return {INSTANCE(Quinlan1999Order8A,
                   0.07 * Second,
                   1.00044972306534419e-13 * Metre,
                   1.00587940754515159e-13 * Metre / Second,
                   3.93946996135596805e-08 * Joule),
          INSTANCE(Quinlan1999Order8B,
                   0.055 * Second,
                   9.97882332320898513e-14 * Metre,
                   1.00953967407946266e-13 * Metre / Second,
                   2.19802622769549316e-08 * Joule),
          INSTANCE(QuinlanTremaine1990Order8,
                   0.3 * Second,
                   9.98298665955132947e-14 * Metre,
                   1.00752739484732956e-13 * Metre / Second,
                   6.42628611435824837e-08 * Joule),
          INSTANCE(QuinlanTremaine1990Order10,
                   0.3 * Second,
                   9.96980276113390573e-14 * Metre,
                   1.02442360150334366e-13 * Metre / Second,
                   1.03418451580239434e-09 * Joule),
          INSTANCE(QuinlanTremaine1990Order12,
                   0.21 * Second,
                   9.90457715843717779e-14 * Metre,
                   1.05165876007617953e-13 * Metre / Second,
                   4.14703826834283973e-11 * Joule)};
}

}  // namespace

class SymmetricLinearMultistepIntegratorTest
    : public ::testing::TestWithParam<SimpleHarmonicMotionTestInstance> {
 public:
  SymmetricLinearMultistepIntegratorTest() {
    google::LogToStderr();
  }
};

INSTANTIATE_TEST_SUITE_P(SymmetricLinearMultistepIntegratorTests,
                        SymmetricLinearMultistepIntegratorTest,
                        ValuesIn(Instances()));

// Tests that the error in energy does not correlate with the number of steps
// taken, and checks the actual value of the error.
TEST_P(SymmetricLinearMultistepIntegratorTest, Symplecticity) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 500 * Second;
  Time const step = 0.2 * Second;

  Mass const m = 1 * Kilogram;
  Stiffness const k = si::Unit<Stiffness>;
  Energy const initial_energy =
      0.5 * m * Pow<2>(v_initial) + 0.5 * k * Pow<2>(q_initial);

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };

  auto const instance =
      GetParam().integrator.NewInstance(problem, append_state, step);
  EXPECT_OK(instance->Solve(t_final));

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
  EXPECT_THAT(correlation, Lt(0.011));
  Power const slope = Slope(time, energy_error);
  LOG(INFO) << "Slope                                     : " << slope;
  EXPECT_THAT(Abs(slope), Lt(2e-6 * si::Unit<Power>));
  LOG(INFO) << "Maximum energy error                      : " <<
      max_energy_error;
  EXPECT_EQ(GetParam().expected_energy_error, max_energy_error);
}

// Integrates with diminishing step sizes, and checks the order of convergence.
TEST_P(SymmetricLinearMultistepIntegratorTest, Convergence) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 100 * Second;

  Time step = GetParam().beginning_of_convergence;
  int const step_sizes = 50;
  double const step_reduction = 1.1;
  std::vector<double> log_step_sizes;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_q_errors;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_p_errors;
  log_step_sizes.reserve(step_sizes);

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  ODE::State final_state;
  auto const append_state = [&final_state](ODE::State const& state) {
    final_state = state;
  };

  Logger logger(
      TEMP_DIR / ("convergence." + GetParam().name + ".generated.wl"),
      /*make_unique=*/false);
  for (int i = 0; i < step_sizes; ++i, step /= step_reduction) {
    auto const instance =
        GetParam().integrator.NewInstance(problem, append_state, step);
    EXPECT_OK(instance->Solve(t_final));
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
    logger.Append("logStepSizes", log_step_sizes.back());
    logger.Append("logQErrors", log_q_errors.back());
    logger.Append("logPErrors", log_p_errors.back());
  }

  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;

#if !defined(_DEBUG)
  EXPECT_THAT(RelativeError(GetParam().order, q_convergence_order),
              Lt(0.02));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.9997), Le(1)));
#endif
  double const v_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_order;
  LOG(INFO) << "Correlation            : " << v_correlation;
#if !defined(_DEBUG)
  EXPECT_THAT(RelativeError(GetParam().order, v_convergence_order), Lt(0.02));
  CHECK_GE(1, v_correlation);
  EXPECT_THAT(v_correlation, AllOf(Gt(0.99993), Le(1)));
#endif
}

TEST_P(SymmetricLinearMultistepIntegratorTest, Termination) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 1630 * Second;
  Time const step = 42 * Second;
  int const steps = static_cast<int>(std::floor((t_final - t_initial) / step));

  int evaluations = 0;

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, &evaluations);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };

  auto const instance =
      GetParam().integrator.NewInstance(problem, append_state, step);
  EXPECT_OK(instance->Solve(t_final));

  EXPECT_EQ(steps, solution.size());
  EXPECT_THAT(solution.back().time.value,
              AllOf(Gt(t_final - step), Le(t_final)));
  for (int i = 0; i < steps; ++i) {
    Time const t = solution[i].time.value - t_initial;
    EXPECT_THAT(t, AlmostEquals((i + 1) * step, 0));
  }
}

// Long integration, change detector.  Also tests the number of steps, their
// spacing, and the number of evaluations.
TEST_P(SymmetricLinearMultistepIntegratorTest, LongIntegration) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
#if defined(_DEBUG)
  Instant const t_final = t_initial + 1 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 1 * Milli(Second);
  int const steps = static_cast<int>((t_final - t_initial) / step) - 1;

  int evaluations = 0;

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, &evaluations);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };

  auto const instance =
      GetParam().integrator.NewInstance(problem, append_state, step);
  EXPECT_OK(instance->Solve(t_final));

  EXPECT_EQ(steps, solution.size());
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
#if !defined(_DEBUG)
  EXPECT_EQ(GetParam().expected_position_error, q_error);
  EXPECT_EQ(GetParam().expected_velocity_error, v_error);
#endif
}

// Tests that serialization and deserialization work.
TEST_P(SymmetricLinearMultistepIntegratorTest, Serialization) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Time const step = 0.2 * Second;

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {t_initial, {q_initial}, {v_initial}};
  auto const append_state = [&solution](ODE::State const& state) {
    solution.push_back(state);
  };

  auto const instance1 =
      GetParam().integrator.NewInstance(problem, append_state, step);
  serialization::IntegratorInstance message1;
  instance1->WriteToMessage(&message1);
  auto const instance2 =
      FixedStepSizeIntegrator<ODE>::Instance::ReadFromMessage(
          message1, harmonic_oscillator, append_state);
  serialization::IntegratorInstance message2;
  instance2->WriteToMessage(&message2);
  EXPECT_THAT(message1, EqualsProto(message2));
}

}  // namespace integrators
}  // namespace principia
