#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"  // NOLINT(whitespace/line_length)

#include <algorithm>
#include <limits>
#include <vector>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/legendre.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/integration.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_generalized_runge_kutta_nyström_integrator {  // NOLINT(whitespace/line_length)

using numerics::EstrinEvaluator;
using numerics::LegendrePolynomial;
using quantities::Abs;
using quantities::Time;
using quantities::Variation;
using quantities::si::Centi;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Newton;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::ComputeЧебышёвPolynomialSecondDerivative;
using testing_utilities::ComputeLegendrePolynomialSecondDerivative;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::operator""_⑴;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::std::placeholders::_4;
using ::testing::ElementsAreArray;
using ::testing::Lt;

using ODE = ExplicitSecondOrderOrdinaryDifferentialEquation<double>;

namespace {

double ToleranceToErrorRatio(
    Time const& h,
    ODE::SystemStateError const& error,
    double const& tolerance,
    Variation<double> const& derivative_tolerance,
    std::function<void(bool tolerable)> callback) {
  double const r =
      std::min(tolerance / Abs(error.position_error[0]),
               derivative_tolerance / Abs(error.velocity_error[0]));
  callback(r > 1.0);
  return r;
}

}  // namespace

class EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegratorTest
    : public ::testing::Test {};

TEST_F(EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegratorTest, Legendre) {
  AdaptiveStepSizeIntegrator<ODE> const& integrator =
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
          methods::Fine1987RKNG34,
          double>();
  // TODO(egg): Change that back to 15 once compiling LegendrePolynomial<15>
  // becomes tractable.
  constexpr int degree = 3;
  double const x_initial = 0;
  Variation<double> const v_initial =
      -3 / (2 * Second);  // -6435 / (2048 * Second);
  Instant const t_initial;
  Instant const t_final = t_initial + 0.99 * Second;
  double const tolerance = 1e-6;
  Variation<double> const derivative_tolerance = 1e-6 / Second;

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
  ODE legendre_equation;
  legendre_equation.compute_acceleration =
      std::bind(ComputeLegendrePolynomialSecondDerivative<degree>,
                _1, _2, _3, _4, &evaluations);
  IntegrationProblem<ODE> problem;
  problem.equation = legendre_equation;
  problem.initial_state = {{x_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.push_back(state);
  };

  AdaptiveStepSizeIntegrator<ODE>::Parameters const parameters(
      /*first_time_step=*/t_final - t_initial,
      /*safety_factor=*/0.9);
  auto const tolerance_to_error_ratio = std::bind(ToleranceToErrorRatio,
                                                  _1,
                                                  _2,
                                                  tolerance,
                                                  derivative_tolerance,
                                                  step_size_callback);
  auto instance = integrator.NewInstance(
      problem, append_state, tolerance_to_error_ratio, parameters);
  auto outcome = instance->Solve(t_final);
  EXPECT_EQ(termination_condition::Done, outcome.error());

  double max_error{};
  Variation<double> max_derivative_error{};
  for (ODE::SystemState const& state : solution) {
    double const x = (state.time.value - t_initial) / (1 * Second);
    double const error =
        AbsoluteError(LegendrePolynomial<degree, EstrinEvaluator>().Evaluate(x),
                      state.positions[0].value);
    Variation<double> const derivative_error = AbsoluteError(
        LegendrePolynomial<degree, EstrinEvaluator>().Derivative().Evaluate(x) /
            (1 * Second),
        state.velocities[0].value);
    max_error = std::max(max_error, error);
    max_derivative_error = std::max(max_derivative_error, derivative_error);
  }
  EXPECT_THAT(max_error, IsNear(172e-6_⑴));
  EXPECT_THAT(max_derivative_error, IsNear(4.54e-3_⑴ / Second));
}

}  // namespace internal_embedded_explicit_generalized_runge_kutta_nyström_integrator  // NOLINT
}  // namespace integrators
}  // namespace principia
