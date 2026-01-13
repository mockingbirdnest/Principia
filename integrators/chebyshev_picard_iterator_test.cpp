#include "integrators/chebyshev_picard_iterator.hpp"

#include "geometry/instant.hpp"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace integrators {

using ::testing::Values;

using namespace principia::geometry::_instant;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_chebyshev_picard_iterator;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

using ODE = ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length>;

namespace {

struct ChebyshevPicardIteratorTestParam {
  // The step size used for this test.
  ODE::IndependentVariableDifference step;

  // The stopping criterion used for this test.
  double stopping_criterion;

  // The number of ULPs the integrated solution is allowed to vary from the
  // analytic solution.
  int ulps;
};

class ChebyshevPicardIteratorTest
    : public testing::TestWithParam<ChebyshevPicardIteratorTestParam> {};

TEST_P(ChebyshevPicardIteratorTest, LinearConvergence) {
  // The first order ODE y′ = y.
  //
  // The solution to this is known; it is y = Ce^t (for y(0) = 1, which we use
  // in this test, C = 1).
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
  Time const step = GetParam().step;
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
  ChebyshevPicardIterator<ODE> const integrator(ChebyshevPicardIterationParams{
      .M = 64,
      .N = 64,
      .max_iterations = 64,
      .stopping_criterion = GetParam().stopping_criterion,
  });

  auto const instance = integrator.NewInstance(problem, append_state, step);
  EXPECT_OK(instance->Solve(t_final));

  // Verify the results are close to the known solution.
  for (const auto& state : states) {
    double t = (state.s.value - Instant()) / Second;
    double y = std::get<0>(state.y).value / Metre;
    EXPECT_THAT(y, AlmostEquals(std::exp(t), 0, GetParam().ulps)) << "t=" << t;
  }
}

// Although in theory the iteration for y′ = y should converge for intervals of
// length < 40 (for sufficiently high N), in practice the convergence degrades
// far earlier.
INSTANTIATE_TEST_SUITE_P(AllChebyshevPicardIteratorTests,
                         ChebyshevPicardIteratorTest,
                         Values(
                             ChebyshevPicardIteratorTestParam{
                                 .stopping_criterion = 1e-16,
                                 .step = 1 * Second,
                                 .ulps = 2,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .stopping_criterion = 1e-16,
                                 .step = 2 * Second,
                                 .ulps = 5,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .stopping_criterion = 1e-16,
                                 .step = 4 * Second,
                                 .ulps = 15,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .stopping_criterion = 1e-15,
                                 .step = 8 * Second,
                                 .ulps = (int)1.2e3,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .stopping_criterion = 1e-11,
                                 .step = 16 * Second,
                                 .ulps = (int)4e7,
                             }));

}  // namespace

}  // namespace integrators
}  // namespace principia