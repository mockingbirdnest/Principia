#include "integrators/explicit_linear_multistep_integrator.hpp"

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "base/macros.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/logger.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_runge_kutta_integrator {

using geometry::Instant;
using quantities::Length;
using quantities::Mass;
using quantities::SpecificImpulse;
using quantities::Speed;
using quantities::Time;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::ApproximateQuantity;
using testing_utilities::IsNear;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;
using testing_utilities::operator""_;
using ::testing::ValuesIn;

// The execution time is exponential in |step_sizes|.
constexpr int step_sizes = 110;
constexpr double step_reduction = 1.1;

#define PARAM(integrator,                                            \
              initial_number_of_steps,                               \
              expected_q_convergence,                                \
              expected_q_correlation,                                \
              expected_v_convergence,                                \
              expected_v_correlation)                                \
  IntegratorTestParam(                                               \
      ExplicitLinearMultistepIntegrator<methods::integrator, ODE>(), \
      #integrator,                                                   \
      (initial_number_of_steps),                                     \
      (expected_q_convergence),                                      \
      (expected_q_correlation),                                      \
      (expected_v_convergence),                                      \
      (expected_v_correlation))

using ODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length, Speed>;

struct IntegratorTestParam final {
  template<typename Integrator>
  IntegratorTestParam(
      Integrator const& integrator,
      std::string const& name,
      int const initial_number_of_steps,
      ApproximateQuantity<double> const& expected_q_convergence,
      ApproximateQuantity<double> const& expected_q_correlation,
      ApproximateQuantity<double> const& expected_v_convergence,
      ApproximateQuantity<double> const& expected_v_correlation)
      : integrator(integrator),
        name(name),
        initial_number_of_steps(initial_number_of_steps),
        expected_q_convergence(expected_q_convergence),
        expected_q_correlation(expected_q_correlation),
        expected_v_convergence(expected_v_convergence),
        expected_v_correlation(expected_v_correlation) {}

  FixedStepSizeIntegrator<ODE> const& integrator;
  std::string const name;
  int const initial_number_of_steps;
  ApproximateQuantity<double> const expected_q_convergence;
  ApproximateQuantity<double> const expected_q_correlation;
  ApproximateQuantity<double> const expected_v_convergence;
  ApproximateQuantity<double> const expected_v_correlation;
};

// This allows the test output to be legible, i.e.,
// "where GetParam() = Leapfrog" rather than
// "where GetParam() = n-byte object <hex>"
std::ostream& operator<<(std::ostream& stream,
                         IntegratorTestParam const& param) {
  return stream << param.name;
}

std::vector<IntegratorTestParam> IntegratorTestParams() {
  // The |initial_number_of_steps| below were carefully chosen using
  // Mathematica to only select the domain where |p| and |step| are properly
  // correlated.
  return {PARAM(AdamsBashforthOrder2,
                10,
                1.964_(1),
                0.99985_(1),
                1.942_(1),
                0.99963_(1)),
          PARAM(AdamsBashforthOrder3,
                20,
                2.928_(1),
                0.99980_(1),
                2.898_(1),
                0.99960_(1)),
          PARAM(AdamsBashforthOrder4,
                40,
                3.86_(1),
                0.99977_(1),
                3.81_(1),
                0.99960_(1)),
          PARAM(AdamsBashforthOrder5,
                70,
                4.78_(1),
                0.99981_(1),
                4.72_(1),
                0.99971_(1))};
}

class ExplicitLinearMultistepIntegratorTest
      : public ::testing::TestWithParam<IntegratorTestParam> {
 public:
  ExplicitLinearMultistepIntegratorTest() {
    google::LogToStderr();
  }
};

#if !defined(_DEBUG)
INSTANTIATE_TEST_SUITE_P(ExplicitLinearMultistepIntegratorTests,
                         ExplicitLinearMultistepIntegratorTest,
                         ValuesIn(IntegratorTestParams()));

// Integrates with diminishing step sizes, and checks the order of convergence.
TEST_P(ExplicitLinearMultistepIntegratorTest, Convergence) {
  LOG(INFO) << GetParam();
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
  // Before the singularity.
  Instant const t_final = t_initial + 0.9 * initial_mass / mass_flow;
  Length const q_initial = 0 * Metre;
  Speed const v_initial = 0 * Metre / Second;
  auto const mass = [initial_mass, t_initial, mass_flow](Instant const& t) {
    return initial_mass - (t - t_initial) * mass_flow;
  };

  Time step = (t_final - t_initial) / GetParam().initial_number_of_steps;
  std::vector<double> log_step_sizes;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_q_errors;
  log_step_sizes.reserve(step_sizes);
  std::vector<double> log_p_errors;
  log_step_sizes.reserve(step_sizes);

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
  ODE::State final_state;
  auto const append_state = [&final_state](ODE::State const& state) {
    final_state = state;
  };
  problem.equation = rocket_equation;
  problem.initial_state = {t_initial, {{q_initial}, {v_initial}}};

  FixedStepSizeIntegrator<ODE> const& integrator = GetParam().integrator;

  mathematica::Logger logger(
      TEMP_DIR / ("convergence." + GetParam().name + ".generated.wl"),
      /*make_unique=*/false);
  for (int i = 0; i < step_sizes; ++i, step /= step_reduction) {
    auto const instance =
        integrator.NewInstance(problem, append_state, step);
    EXPECT_OK(instance->Solve(t_final));
    Time const t = final_state.s.value - t_initial;
    Length const& q = std::get<0>(final_state.y).value;
    Speed const& v = std::get<1>(final_state.y).value;
    double const log_q_error = std::log10(RelativeError(
        q,
        specific_impulse *
            (t + (t - initial_mass / mass_flow) *
                     std::log(initial_mass / mass(final_state.s.value)))));
    double const log_p_error = std::log10(RelativeError(
        v,
        specific_impulse * std::log(initial_mass / mass(final_state.s.value))));
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
  double const q_convergence_error = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_error;
  LOG(INFO) << "Correlation            : " << q_correlation;

  EXPECT_THAT(q_convergence_error, IsNear(GetParam().expected_q_convergence));
  EXPECT_THAT(q_correlation, IsNear(GetParam().expected_q_correlation));
  double const v_convergence_error = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_error;
  LOG(INFO) << "Correlation            : " << v_correlation;
  EXPECT_THAT(v_convergence_error, IsNear(GetParam().expected_v_convergence));
  EXPECT_THAT(v_correlation, IsNear(GetParam().expected_v_correlation));
}
#endif

}  // namespace internal_explicit_runge_kutta_integrator
}  // namespace integrators
}  // namespace principia
