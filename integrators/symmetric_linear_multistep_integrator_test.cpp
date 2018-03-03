
#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include <algorithm>
#include <vector>
#include <string>

#include "base/file.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using base::OFStream;
using geometry::Instant;
using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Energy;
using quantities::Length;
using quantities::Mass;
using quantities::Pow;
using quantities::Power;
using quantities::Sin;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Stiffness;
using quantities::Time;
using quantities::si::Joule;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::ComputeHarmonicOscillatorAcceleration1D;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
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

template<typename Integrator>
void TestTermination(Integrator const& integrator) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 1630 * Second;
  Time const step = 42 * Second;
  int const steps = static_cast<int>(std::floor((t_final - t_initial) / step));

  int evaluations = 0;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);
  instance->Solve(t_final);

  EXPECT_EQ(steps, solution.size());
  EXPECT_THAT(solution.back().time.value,
              AllOf(Gt(t_final - step), Le(t_final)));
  Length q_error;
  Speed v_error;
  for (int i = 0; i < steps; ++i) {
    Time const t = solution[i].time.value - t_initial;
    EXPECT_THAT(t, AlmostEquals((i + 1) * step, 0));
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
#if defined(_DEBUG)
  Instant const t_final = t_initial + 1 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 1 * Milli(Second);
  int const steps = static_cast<int>((t_final - t_initial) / step) - 1;

  int evaluations = 0;

  std::vector<ODE::SystemState> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);
  instance->Solve(t_final);

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
  EXPECT_EQ(expected_position_error, q_error);
  EXPECT_EQ(expected_velocity_error, v_error);
#endif
}

// Integrates with diminishing step sizes, and checks the order of convergence.
template<typename Integrator>
void TestConvergence(Integrator const& integrator,
                     std::string const& name,
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
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{q_initial}, {v_initial}, t_initial};
  ODE::SystemState final_state;
  auto const append_state = [&final_state](ODE::SystemState const& state) {
    final_state = state;
  };

  for (int i = 0; i < step_sizes; ++i, step /= step_reduction) {
    auto const instance = integrator.NewInstance(problem, append_state, step);
    instance->Solve(t_final);
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

  {
    std::experimental::filesystem::path filename;
    filename += "convergence.";
    filename += name;
    filename += ".generated.wl";
    OFStream file(TEMP_DIR / filename);
    file << mathematica::Assign("logStepSizes", log_step_sizes);
    file << mathematica::Assign("logQErrors", log_q_errors);
    file << mathematica::Assign("logPErrors", log_p_errors);
  }

  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;

#if !defined(_DEBUG)
  EXPECT_THAT(RelativeError(integrator.order, q_convergence_order),
              Lt(0.02));
  EXPECT_THAT(q_correlation, IsNear(1.0, /*tolerance=*/1.0005));
#endif
  double const v_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_order;
  LOG(INFO) << "Correlation            : " << v_correlation;
#if !defined(_DEBUG)
  EXPECT_THAT(RelativeError(integrator.order, v_convergence_order), Lt(0.02));
  EXPECT_THAT(v_correlation, IsNear(1.0, /*tolerance=*/1.0002));
#endif
}

// Tests that the error in energy does not correlate with the number of steps
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
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);
  instance->Solve(t_final);

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
  EXPECT_THAT(Abs(slope), Lt(2e-6 * SIUnit<Power>()));
  LOG(INFO) << "Maximum energy error                      : " <<
      max_energy_error;
  EXPECT_EQ(expected_energy_error, max_energy_error);
}

// Tests that serialization and deserialization work.
template<typename Integrator>
void TestSerialization(Integrator const& integrator) {
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Time const step = 0.2 * Second;

  Mass const m = 1 * Kilogram;
  Stiffness const k = SIUnit<Stiffness>();
  Energy const initial_energy =
      0.5 * m * Pow<2>(v_initial) + 0.5 * k * Pow<2>(q_initial);

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

  auto const instance1 = integrator.NewInstance(problem, append_state, step);
  serialization::IntegratorInstance message1;
  instance1->WriteToMessage(&message1);
  auto const instance2 = FixedStepSizeIntegrator<
      typename Integrator::ODE>::Instance::ReadFromMessage(message1,
                                                           harmonic_oscillator,
                                                           append_state);
  serialization::IntegratorInstance message2;
  instance2->WriteToMessage(&message2);
  EXPECT_THAT(message1, EqualsProto(message2));
}

class SimpleHarmonicMotionTestInstance final {
 public:
  template<typename Integrator>
  SimpleHarmonicMotionTestInstance(Integrator const& integrator,
                                   std::string const& name,
                                   Time const& beginning_of_convergence,
                                   Length const& expected_position_error,
                                   Speed const& expected_velocity_error,
                                   Energy const& expected_energy_error)
      : test_termination_(std::bind(TestTermination<Integrator>, integrator)),
        test_1000_seconds_at_1_millisecond_(
            std::bind(Test1000SecondsAt1Millisecond<Integrator>,
                      integrator,
                      expected_position_error,
                      expected_velocity_error)),
        test_convergence_(
            std::bind(TestConvergence<Integrator>,
                      integrator,
                      name,
                      beginning_of_convergence)),
        test_symplecticity_(
            std::bind(TestSymplecticity<Integrator>,
                      integrator,
                      expected_energy_error)),
        test_serialization_(
            std::bind(TestSerialization<Integrator>,
                      integrator)),
        name_(name) {}

  std::string const& name() const {
    return name_;
  }

  void RunTermination() const {
    test_termination_();
  }

  void Run1000SecondsAt1Millisecond() const {
    test_1000_seconds_at_1_millisecond_();
  }

  void RunConvergence() const {
    test_convergence_();
  }

  void RunSymplecticity() const {
    test_symplecticity_();
  }

  void RunSerialization() const {
    test_serialization_();
  }

 private:
  std::function<void()> test_termination_;
  std::function<void()> test_1000_seconds_at_1_millisecond_;
  std::function<void()> test_convergence_;
  std::function<void()> test_symplecticity_;
  std::function<void()> test_serialization_;
  std::string name_;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         SimpleHarmonicMotionTestInstance instance) {
  return stream << instance.name();
}

// Not testing QuinlanTremaine1990Order14 as its characteristics cannot even be
// computed.
std::vector<SimpleHarmonicMotionTestInstance> Instances() {
  // The |beginning_of_convergence| below were carefully chosen using
  // Mathematica to only select the domain where |p| and |step| are properly
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

INSTANTIATE_TEST_CASE_P(SymmetricLinearMultistepIntegratorTests,
                        SymmetricLinearMultistepIntegratorTest,
                        ValuesIn(Instances()));

TEST_P(SymmetricLinearMultistepIntegratorTest, Symplecticity) {
  LOG(INFO) << GetParam();
  GetParam().RunSymplecticity();
}

TEST_P(SymmetricLinearMultistepIntegratorTest, Convergence) {
  LOG(INFO) << GetParam();
  GetParam().RunConvergence();
}

TEST_P(SymmetricLinearMultistepIntegratorTest, Termination) {
  LOG(INFO) << GetParam();
  GetParam().RunTermination();
}

TEST_P(SymmetricLinearMultistepIntegratorTest, LongIntegration) {
  LOG(INFO) << GetParam();
  GetParam().Run1000SecondsAt1Millisecond();
}

TEST_P(SymmetricLinearMultistepIntegratorTest, Serialization) {
  LOG(INFO) << GetParam();
  GetParam().RunSerialization();
}

}  // namespace integrators
}  // namespace principia
