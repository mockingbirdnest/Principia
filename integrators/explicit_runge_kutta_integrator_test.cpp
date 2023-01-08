#include "integrators/explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "base/macros.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
using testing_utilities::AbsoluteError;
using testing_utilities::IsNear;
using testing_utilities::PearsonProductMomentCorrelationCoefficient;
using testing_utilities::RelativeError;
using testing_utilities::Slope;
using testing_utilities::operator""_;

using ODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length, Speed>;

class ExplicitRungeKuttaIntegratorTest : public ::testing::Test {
 public:
  ExplicitRungeKuttaIntegratorTest() {
    google::LogToStderr();
  }
};

// Integrates with diminishing step sizes, and checks the order of convergence.
TEST_F(ExplicitRungeKuttaIntegratorTest, Convergence) {
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

  Time step = (t_final - t_initial) / 10;
  int const step_sizes = 50;
  double const step_reduction = 1.1;
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
    qʹ[0] = v[0];
    vʹ[0] = mass_flow * specific_impulse / mass(t);
    return absl::OkStatus();
  };
  InitialValueProblem<ODE> problem;
  ODE::State final_state;
  auto const append_state = [&final_state](ODE::State const& state) {
    final_state = state;
  };
  problem.equation = rocket_equation;
  problem.initial_state = {t_initial, {{q_initial}, {v_initial}}};

  FixedStepSizeIntegrator<ODE> const& integrator =
      ExplicitRungeKuttaIntegrator<methods::Kutta1901Vσ1, ODE>();

  for (int i = 0; i < step_sizes; ++i, step /= step_reduction) {
    auto const instance =
        integrator.NewInstance(problem, append_state, step);
    EXPECT_OK(instance->Solve(t_final));
    Time const t = final_state.s.value - t_initial;
    Length const& q = std::get<0>(final_state.y)[0].value;
    Speed const& v = std::get<1>(final_state.y)[0].value;
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
  }
  double const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  double const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;

#if !defined(_DEBUG)
  EXPECT_THAT(AbsoluteError(static_cast<double>(methods::Kutta1901Vσ1::order),
                            q_convergence_order),
              IsNear(0.15_(1)));
  EXPECT_THAT(q_correlation, IsNear(0.9996_(1)));
#endif
  double const v_convergence_order = Slope(log_step_sizes, log_p_errors);
  double const v_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << v_convergence_order;
  LOG(INFO) << "Correlation            : " << v_correlation;
#if !defined(_DEBUG)
  EXPECT_THAT(AbsoluteError(methods::Kutta1901Vσ1::order, v_convergence_order),
              IsNear(0.19_(1)));
  EXPECT_THAT(v_correlation, IsNear(0.9992_(1)));
#endif
}

}  // namespace internal_ordinary_differential_equations

}  // namespace integrators
}  // namespace principia
