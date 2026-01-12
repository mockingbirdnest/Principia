#include "integrators/chebyshev_picard_iterator.hpp"

#include "geometry/instant.hpp"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace integrators {

using ::testing::DoubleNear;

using namespace principia::geometry::_instant;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_chebyshev_picard_iterator;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

using ODE = ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length>;

namespace {

TEST(ChebyshevPicardIteratorTest, LinearConvergence) {
  // The first order ODE y′ = y.
  //
  // The solution to this is known; it is y = Ce^t (for y(0) = 1, which we use
  // in this test, C = 1).
  // The convergence interval (of the Picard iteration) is also known: for order
  // > 40, the convergence interval length is also ~40 (see Bai's dissertation,
  // equation 3.12).
  ODE linear_ode;
  linear_ode.compute_derivative =
      [](Instant const& t, ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = y / Second;
        return absl::OkStatus();
      };

  // Set up the initial value problem.
  Time const step = 1 * Second;
  Instant const t_initial;
  Instant const t_final = t_initial + step;

  InitialValueProblem<ODE> problem;
  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  problem.equation = linear_ode;
  problem.initial_state = {t_initial, {1 * Metre}};

  // Build the integrator and solve the problem.
  const int order = 64;
  ChebyshevPicardIterator<ODE> const integrator(/*order=*/order,
                                                /*sample_points=*/order,
                                                /*max_iterations=*/32,
                                                /*stopping_criterion=*/1e-16);

  auto const instance = integrator.NewInstance(problem, append_state, step);
  EXPECT_OK(instance->Solve(t_final));

  // Verify the results are close to the known solution.
  for (const auto& state : states) {
    double t = (state.s.value - Instant()) / Second;
    double y = std::get<0>(state.y).value / Metre;
    EXPECT_THAT(y, DoubleNear(std::exp(t), 1e-15)) << "t=" << t;
  }
}

}  // namespace

}  // namespace integrators
}  // namespace principia