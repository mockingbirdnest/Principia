#include "integrators/chebyshev_picard_iterator.hpp"

#include "geometry/instant.hpp"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace integrators {

using ::testing::TestWithParam;
using ::testing::Values;

using namespace principia::geometry::_instant;
using namespace principia::integrators::_chebyshev_picard_iterator;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

using ODE = ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Length>;

namespace {

// An initial value problem with a known solution.
struct SolvedInitialValueProblem {
  ODE::IndependentVariable t₀() const { return problem.initial_state.s.value; }

  InitialValueProblem<ODE> problem;
  std::function<ODE::DependentVariables(ODE::IndependentVariable const&)>
      solution;
};

// The first order ODE y′ = y with y₀ = 1.
//
// The solution is y = eᵗ.
SolvedInitialValueProblem LinearProblem() {
  ODE linear_ode;
  linear_ode.compute_derivative =
      [](Instant const& t, ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = y / Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = linear_ode;
  problem.initial_state = {Instant(), {1 * Metre}};

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [](Instant const& t) -> std::tuple<Length> {
        return std::exp((t - Instant()) / Second) * Metre;
      }};
}

// The first order ODE y′ = π/8 (1 + y²) with y(-1) = 0.
//
// The solution is y = tan(π/8 (t + 1)).
SolvedInitialValueProblem TangentProblem() {
  ODE ode;
  ode.compute_derivative =
      [](Instant const& t, ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = π / 8 * (1 + y * y / (Metre * Metre)) * Metre / Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  // TODO(rnlahaye): There is some issue when we start at t=-1. Fix it.
  problem.initial_state = {Instant() - 0.9 * Second,
                           {Tan(π / 8 * 0.1 * Radian) * Metre}};

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [](Instant const& t) -> std::tuple<Length> {
        return Tan(π / 8 * ((t - Instant()) / Second + 1) * Radian) * Metre;
      }};
}

struct ChebyshevPicardIteratorTestParam {
  // The ODE (with solution) under test.
  SolvedInitialValueProblem problem;

  // The Chebyshev polynomial order.
  int N;

  // The step size used for this test.
  ODE::IndependentVariableDifference step;

  // The stopping criterion used for this test.
  double stopping_criterion;

  // The number of ULPs the integrated solution is allowed to vary from the
  // analytic solution.
  int ulps;
};

class ChebyshevPicardIteratorTest
    : public TestWithParam<ChebyshevPicardIteratorTestParam> {};

TEST_P(ChebyshevPicardIteratorTest, LinearConvergence) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = GetParam().problem;
  Time const step = GetParam().step;

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ChebyshevPicardIterator<ODE> const integrator(ChebyshevPicardIterationParams{
      .M = GetParam().N,
      .N = GetParam().N,
      .max_iterations = 64,
      .stopping_criterion = GetParam().stopping_criterion,
  });

  auto const instance =
      integrator.NewInstance(problem.problem, append_state, step);
  EXPECT_OK(instance->Solve(problem.t₀() + step));

  // Verify the results are close to the known solution.
  for (const auto& state : states) {
    auto t = state.s.value;
    auto y = std::get<0>(state.y).value;
    EXPECT_THAT(
        y, AlmostEquals(std::get<0>(problem.solution(t)), 0, GetParam().ulps))
        << "t=" << t;
  }
}

// Although in theory the iteration for y′ = y should converge for
// intervals
// of length < 40 for sufficiently high N (see Bai's thesis), in
// practice the convergence degrades far earlier.
INSTANTIATE_TEST_SUITE_P(Linear, ChebyshevPicardIteratorTest,
                         Values(
                             ChebyshevPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-16,
                                 .step = 1 * Second,
                                 .ulps = 2,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-16,
                                 .step = 2 * Second,
                                 .ulps = 5,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-16,
                                 .step = 4 * Second,
                                 .ulps = 15,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-15,
                                 .step = 8 * Second,
                                 .ulps = (int)1.2e3,
                             },
                             ChebyshevPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-11,
                                 .step = 16 * Second,
                                 .ulps = (int)4e7,
                             }));

// This problem appears in the NASA paper.
INSTANTIATE_TEST_SUITE_P(Tangent, ChebyshevPicardIteratorTest,
                         Values(
                             // These are the parameters from the NASA paper.
                             ChebyshevPicardIteratorTestParam{
                                 .problem = TangentProblem(),
                                 .N = 16,
                                 .stopping_criterion = 0.5e-10,
                                 .step = 2 * Second,
                                 .ulps = (int)3.4e5,
                             },
                             // We can achieve better accuracy with higher N and
                             // a more stringent stopping criterion.
                             ChebyshevPicardIteratorTestParam{
                                 .problem = TangentProblem(),
                                 .N = 32,
                                 .stopping_criterion = 1e-16,
                                 .step = 2 * Second,
                                 .ulps = 21,
                             }));

}  // namespace

}  // namespace integrators
}  // namespace principia