#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <vector>
#include <string>

#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/integrators.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
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
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
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

#define INSTANCE(integrator,                                             \
                 beginning_of_convergence,                               \
                 expected_position_error,                                \
                 expected_velocity_error,                                \
                 expected_energy_error)                                  \
  SimpleHarmonicMotionTestInstance(                                      \
      SymplecticRungeKuttaNyströmIntegrator<methods::integrator, ODE>(), \
      #integrator,                                                       \
      (beginning_of_convergence),                                        \
      (expected_position_error),                                         \
      (expected_velocity_error),                                         \
      (expected_energy_error),                                           \
      true)

#define SPRK_INSTANCE(integrator,                                         \
                      composition,                                        \
                      beginning_of_convergence,                           \
                      expected_position_error,                            \
                      expected_velocity_error,                            \
                      expected_energy_error)                              \
  SimpleHarmonicMotionTestInstance(                                       \
      SymplecticRungeKuttaNyströmIntegrator<                              \
          methods::integrator,                                            \
          serialization::FixedStepSizeIntegrator::composition,            \
          ODE>(),                                                         \
      "AsSymplecticRungeKuttaNyström<" #integrator ", " #composition ">", \
      (beginning_of_convergence),                                         \
      (expected_position_error),                                          \
      (expected_velocity_error),                                          \
      (expected_energy_error),                                            \
      false)

using ODE = SpecialSecondOrderDifferentialEquation<Length>;

namespace {

struct SimpleHarmonicMotionTestInstance final {
  template<typename Integrator>
  SimpleHarmonicMotionTestInstance(Integrator const& integrator,
                                   std::string const& name,
                                   Time const& beginning_of_convergence,
                                   Length const& expected_position_error,
                                   Speed const& expected_velocity_error,
                                   Energy const& expected_energy_error,
                                   bool const serializable)
      : integrator(integrator),
        composition(integrator.composition),
        evaluations(integrator.evaluations),
        order(integrator.order),
        time_reversible(integrator.time_reversible),
        name(name),
        beginning_of_convergence(beginning_of_convergence),
        expected_position_error(expected_position_error),
        expected_velocity_error(expected_velocity_error),
        expected_energy_error(expected_energy_error),
        serializable(serializable) {}

  FixedStepSizeIntegrator<ODE> const& integrator;
  serialization::FixedStepSizeIntegrator::CompositionMethod const composition;
  int const evaluations;
  int const order;
  bool const time_reversible;

  std::string const name;

  Time const beginning_of_convergence;
  Length const expected_position_error;
  Speed const expected_velocity_error;
  Energy const expected_energy_error;

  bool const serializable;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         SimpleHarmonicMotionTestInstance instance) {
  return stream << instance.name;
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
                   +6.29052365752613700e-13 * Joule),
          // We test `NewtonDelambreStørmerVerletLeapfrog` both as `ABA` and
          // `BAB` (sometimes called leapfrog and pseudo-leapfrog) for coverage.
          // We test the others as BAB integrators only.
          SPRK_INSTANCE(NewtonDelambreStørmerVerletLeapfrog,
                        ABA,
                        0.4 * Second,
                        +4.15606749774469295e-05 * Metre,
                        +4.16264386218978197e-05 * Metre / Second,
                        +5.05049535215751355e-03 * Joule),
          SPRK_INSTANCE(NewtonDelambreStørmerVerletLeapfrog,
                        BAB,
                        0.4 * Second,
                        +4.15606749774714325e-05 * Metre,
                        +4.16261832565070940e-05 * Metre / Second,
                        +4.99999039863668893e-03 * Joule),
          SPRK_INSTANCE(Ruth1983,
                        BA,
                        0.1 * Second,
                        +2.77767866216707682e-11 * Metre,
                        +7.01570745942348140e-13 * Metre / Second,
                        +1.15535032619074052e-04 * Joule),
          SPRK_INSTANCE(鈴木1990,
                        BAB,
                        1 * Second,
                        +1.17211795824800902e-12 * Metre,
                        +1.17471483929154630e-12 * Metre / Second,
                        +5.75983521433620638e-06 * Joule),
          SPRK_INSTANCE(吉田1990Order6A,
                        BAB,
                        1.5 * Second,
                        +8.30863156053851526e-14 * Metre,
                        +8.30672336471494077e-14 * Metre / Second,
                        1.28253665132582739e-07 * Joule),
          SPRK_INSTANCE(吉田1990Order6B,
                        BAB,
                        1 * Second,
                        +3.32525673663042198e-13 * Metre,
                        +3.32810168313102395e-13 * Metre / Second,
                        +3.39431978840787352e-06 * Joule),
          SPRK_INSTANCE(吉田1990Order6C,
                        BAB,
                        1 * Second,
                        +9.58053081312471022e-14 * Metre,
                        +9.58625540059543368e-14 * Metre / Second,
                        +3.58056353333413568e-06 * Joule),
          SPRK_INSTANCE(吉田1990Order8A,
                        BAB,
                        0.043 * Second,
                        +6.06330957664269476e-13 * Metre,
                        +6.06924233093053545e-13 * Metre / Second,
                        +1.49030436397135091e-05 * Joule),
          SPRK_INSTANCE(吉田1990Order8B,
                        BAB,
                        0.5 * Second,
                        +4.91648388667442759e-13 * Metre,
                        +4.92134111240716265e-13 * Metre / Second,
                        +1.33083068010186878e-07 * Joule),
          SPRK_INSTANCE(吉田1990Order8C,
                        BAB,
                        0.9 * Second,
                        +3.13037790133918747e-13 * Metre,
                        +3.13291059761411361e-13 * Metre / Second,
                        +4.68151000188044009e-08 * Joule),
          SPRK_INSTANCE(吉田1990Order8D,
                        BAB,
                        1.1 * Second,
                        +2.20309881449054501e-13 * Metre,
                        +2.20490292690556089e-13 * Metre / Second,
                        +1.58094315416690279e-10 * Joule),
          SPRK_INSTANCE(吉田1990Order8E,
                        BAB,
                        0.3 * Second,
                        +1.39072781069060625e-13 * Metre,
                        +1.39159517242859465e-13 * Metre / Second,
                        +3.42872182312881080e-08 * Joule),
          SPRK_INSTANCE(CandyRozmus1991ForestRuth1990,
                        BAB,
                        0.5 * Second,
                        +6.63488881891272086e-11 * Metre,
                        +6.64553828633174248e-11 * Metre / Second,
                        +6.26859072366814374e-05 * Joule),
          SPRK_INSTANCE(McLachlanAtela1992Order2Optimal,
                        BA,
                        0.7 * Second,
                        +2.01685999379921758e-05 * Metre,
                        +2.02003819904818379e-05 * Metre / Second,
                        +8.63068191495619530e-05 * Joule),
          SPRK_INSTANCE(McLachlanAtela1992Order3Optimal,
                        BA,
                        0.09 * Second,
                        +1.21425465168800706e-11 * Metre,
                        +3.68977418063742846e-13 * Metre / Second,
                        +4.72513762963533424e-05 * Joule),
          SPRK_INSTANCE(McLachlan1995S2,
                        BAB,
                        0.1 * Second,
                        +1.20001294944783662e-05 * Metre,
                        +1.20190260074261876e-05 * Metre / Second,
                        +4.47986262497312993e-05 * Joule),
          SPRK_INSTANCE(McLachlan1995SS5,
                        BAB,
                        0.1 * Second,
                        +1.99751326590558165e-12 * Metre,
                        +2.00208426226478053e-12 * Metre / Second,
                        +3.99026027320115162e-06 * Joule),
          SPRK_INSTANCE(McLachlan1995S4,
                        BAB,
                        0.1 * Second,
                        +2.20774787340616285e-13 * Metre,
                        +2.21161630675759113e-13 * Metre / Second,
                        +2.11567735630691089e-06 * Joule),
          SPRK_INSTANCE(McLachlan1995S5,
                        BAB,
                        1 * Second,
                        +7.15746906188030607e-14 * Metre,
                        +7.15608128309952463e-14 * Metre / Second,
                        +7.82877597083064813e-07 * Joule),
          SPRK_INSTANCE(McLachlan1995SS9,
                        BAB,
                        1 * Second,
                        +1.23734356094473696e-13 * Metre,
                        +1.23817622821320583e-13 * Metre / Second,
                        +1.29530730030857910e-08 * Joule),
          SPRK_INSTANCE(McLachlan1995SS15,
                        BAB,
                        1 * Second,
                        +9.67281810204667636e-14 * Metre,
                        +9.67836921716980214e-14 * Metre / Second,
                        +1.21325172131037107e-11 * Joule),
          SPRK_INSTANCE(McLachlan1995SS17,
                        BAB,
                        1 * Second,
                        +8.34124436188687923e-14 * Metre,
                        +8.34471380883883285e-14 * Metre / Second,
                        +2.24043006369356590e-12 * Joule),
          SPRK_INSTANCE(BlanesMoan2002S6,
                        BAB,
                        1 * Second,
                        +1.23803745033512769e-13 * Metre,
                        +1.23966809040254589e-13 * Metre / Second,
                        +2.24300146345335349e-07 * Joule),
          SPRK_INSTANCE(BlanesMoan2002S10,
                        BAB,
                        1 * Second,
                        +7.67632485354496907e-14 * Metre,
                        +7.67372276833100386e-14 * Metre / Second,
                        +2.45151621225403460e-10 * Joule)};
}

}  // namespace

class SymplecticRungeKuttaNyströmIntegratorTest
    : public ::testing::TestWithParam<SimpleHarmonicMotionTestInstance> {
 public:
  SymplecticRungeKuttaNyströmIntegratorTest() {
    google::LogToStderr();
  }
};

INSTANTIATE_TEST_SUITE_P(SymplecticRungeKuttaNyströmIntegratorTests,
                        SymplecticRungeKuttaNyströmIntegratorTest,
                        ValuesIn(Instances()));

// If the integrator is time-reversible, checks that integrating back and
// forth with a large time step yields the initial value.  If it is not, checks
// that there is an error when integrating back and forth.
TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, TimeReversibility) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 100 * Second;
  // While time-reversibility is an exact property independent of the equation,
  // if the computed solution explodes, all bits are lost and time-reversibility
  // utterly fails in finite precision.  A step of 1 second would break 吉田's
  // 8A method.
  Time const step = 0.5 * Second;

  std::vector<ODE::State> solution;
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D,
                _1, _2, _3, /*evaluations=*/nullptr);
  InitialValueProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  ODE::State final_state;
  auto const append_state = [&final_state](ODE::State const& state) {
    final_state = state;
  };

  {
    problem.initial_state = {t_initial, {q_initial}, {v_initial}};
    auto const instance =
        GetParam().integrator.NewInstance(problem, append_state, step);
    EXPECT_OK(instance->Solve(t_final));
  }

  {
    problem.initial_state = final_state;
    auto const instance =
        GetParam().integrator.NewInstance(problem, append_state, -step);
    EXPECT_OK(instance->Solve(t_initial));
  }

  EXPECT_EQ(t_initial, final_state.time.value);
  if (GetParam().time_reversible) {
    EXPECT_THAT(final_state.positions[0].value,
                AlmostEquals(q_initial, 0, 30));
    EXPECT_THAT(final_state.velocities[0].value,
                VanishesBefore(v_amplitude, 0, 70));
  } else {
    EXPECT_THAT(AbsoluteError(q_initial,
                              final_state.positions[0].value),
                Gt(1e-6 * Metre));
    EXPECT_THAT(AbsoluteError(v_initial,
                              final_state.velocities[0].value),
                Gt(1e-6 * Metre / Second));
  }
}

// Test that the error in energy does not correlate with the number of steps
// taken, and checks the actual value of the error.
TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Symplecticity) {
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
  EXPECT_THAT(correlation, Lt(2e-3));
  Power const slope = Slope(time, energy_error);
  LOG(INFO) << "Slope                                     : " << slope;
  EXPECT_THAT(Abs(slope), Lt(2e-6 * si::Unit<Power>));
  LOG(INFO) << "Maximum energy error                      : " <<
      max_energy_error;
  EXPECT_EQ(GetParam().expected_energy_error, max_energy_error);
}

// Integrates with diminishing step sizes, and checks the order of convergence.
TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Convergence) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Speed const v_amplitude = 1 * Metre / Second;
  AngularFrequency const ω = 1 * Radian / Second;
  Instant const t_initial;
#if defined(_DEBUG)
  Instant const t_final = t_initial + 10 * Second;
#else
  Instant const t_final = t_initial + 100 * Second;
#endif

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
  }
  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;

#if !defined(_DEBUG)
  EXPECT_THAT(AbsoluteError(static_cast<double>(GetParam().order),
                            q_convergence_order),
              Lt(0.10));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.993), Le(1)));
#endif
  double const v_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_order;
  LOG(INFO) << "Correlation            : " << v_correlation;
#if !defined(_DEBUG)
  // SPRKs with odd convergence order have a higher convergence order in p.
  EXPECT_THAT(AbsoluteError(GetParam().order + (GetParam().order % 2),
                            v_convergence_order),
              Lt(0.16));
  EXPECT_THAT(v_correlation, AllOf(Gt(0.994), Le(1)));
#endif
}

TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Termination) {
  LOG(INFO) << GetParam();
  Length const q_initial = 1 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  Instant const t_initial;
  Instant const t_final = t_initial + 163 * Second;
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
  switch (GetParam().composition) {
    case serialization::FixedStepSizeIntegrator::BA:
    case serialization::FixedStepSizeIntegrator::ABA:
      EXPECT_EQ(steps * GetParam().evaluations, evaluations);
      break;
    case serialization::FixedStepSizeIntegrator::BAB:
      EXPECT_EQ(steps * GetParam().evaluations + 1, evaluations);
      break;
    default:
      LOG(FATAL) << "Invalid composition";
  }
  for (int i = 0; i < steps; ++i) {
    Time const t = solution[i].time.value - t_initial;
    EXPECT_THAT(t, AlmostEquals((i + 1) * step, 0));
  }
}

// Long integration, change detector.  Also tests the number of steps, their
// spacing, and the number of evaluations.
TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, LongIntegration) {
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
  switch (GetParam().composition) {
    case serialization::FixedStepSizeIntegrator::BA:
    case serialization::FixedStepSizeIntegrator::ABA:
      EXPECT_EQ(steps * GetParam().evaluations, evaluations);
      break;
    case serialization::FixedStepSizeIntegrator::BAB:
      EXPECT_EQ(steps * GetParam().evaluations + 1, evaluations);
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
#if !defined(_DEBUG)
  EXPECT_EQ(GetParam().expected_position_error, q_error);
  EXPECT_EQ(GetParam().expected_velocity_error, v_error);
#endif
}

// Tests that serialization and deserialization work.
TEST_P(SymplecticRungeKuttaNyströmIntegratorTest, Serialization) {
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
