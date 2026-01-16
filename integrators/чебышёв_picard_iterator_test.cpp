#include "integrators/чебышёв_picard_iterator.hpp"

#include <tuple>
#include <vector>

#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace integrators {

using ::testing::DoubleNear;
using ::testing::TestWithParam;
using ::testing::Values;

using namespace principia::geometry::_instant;
using namespace principia::integrators::_чебышёв_picard_iterator;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_matchers;

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
  problem.initial_state = {Instant() - 1 * Second, {0 * Metre}};

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [](Instant const& t) -> std::tuple<Length> {
        return Tan(π / 8 * ((t - Instant()) / Second + 1) * Radian) * Metre;
      }};
}

// The first order ODE y′ = cos(t + εy).
//
// The solution ([Fuk97], equations 31 and 32) is
// y = -γt + 2/ε tan⁻¹((β + σ cos φ)/(1 + σ sin φ))
// where
// σ = α(sin φ + β cos φ)
// φ = ½ (1 - γε) t
// α = 2ε/(1 - ε + √(1 - ε²))
// β = 1/(1 + α) tan(εy₀/2)
// γ = ε/(1 + √(1 - ε²))
SolvedInitialValueProblem PerturbedSinusoidProblem(double ε, double y₀) {
  ODE ode;
  ode.compute_derivative =
      [ε](Instant const& t, ODE::DependentVariables const& dependent_variables,
          ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = Cos(((t - Instant()) / Second + ε * y / Metre) * Radian) * Metre /
             Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  problem.initial_state = {Instant(), {y₀ * Metre}};

  const double one_plus_sqrt_one_minus_ε² = 1 + Sqrt(1 - ε * ε);
  const double α = 2 * ε / (one_plus_sqrt_one_minus_ε² - ε);
  const double β = 1 / (1 + α) * Tan(ε * y₀ / 2 * Radian);
  const double γ = ε / one_plus_sqrt_one_minus_ε²;

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [ε, α, β, γ](Instant const& t) -> std::tuple<Length> {
        const Angle φ = 0.5 * (1 - γ * ε) * (t - Instant()) / Second * Radian;
        const double σ = α * (Sin(φ) + β * Cos(φ));

        return (-γ * (t - Instant()) / Second +
                2 / ε * ArcTan((β + σ * Cos(φ)) / (1 + σ * Sin(φ))) / Radian) *
               Metre;
      }};
}

struct ЧебышёвPicardIteratorTestParam {
  // The ODE (with solution) under test.
  SolvedInitialValueProblem problem;

  // The Чебышёв polynomial order.
  int N;

  // The step size used for this test.
  ODE::IndependentVariableDifference step;

  // The stopping criterion used for this test.
  double stopping_criterion;

  // The distance the integrated solution is allowed to vary from the
  // analytic solution.
  double tolerance;
};

class ЧебышёвPicardIteratorTest
    : public TestWithParam<ЧебышёвPicardIteratorTestParam> {};

TEST_P(ЧебышёвPicardIteratorTest, Convergence) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = GetParam().problem;
  Time const step = GetParam().step;

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIterator<ODE> const integrator(ЧебышёвPicardIterationParams{
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
    EXPECT_THAT(y / Metre, DoubleNear(std::get<0>(problem.solution(t)) / Metre,
                                      GetParam().tolerance))
        << "t=" << t;
  }
}

TEST_F(ЧебышёвPicardIteratorTest, MultipleSteps) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = LinearProblem();
  Time const step = 1 * Second;

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIterator<ODE> const integrator(ЧебышёвPicardIterationParams{
      .M = 64,
      .N = 64,
      .max_iterations = 64,
      .stopping_criterion = 1e-16,
  });

  auto const instance =
      integrator.NewInstance(problem.problem, append_state, step);
  auto const t_final = problem.t₀() + 2.5 * step;

  EXPECT_OK(instance->Solve(t_final));

  // Verify we have reached the desired final time.
  EXPECT_GE(states.back().s.value, t_final);

  // Verify the results are close to the known solution.
  for (const auto& state : states) {
    auto t = state.s.value;
    auto y = std::get<0>(state.y).value;
    EXPECT_THAT(y, AlmostEquals(std::get<0>(problem.solution(t)), 0, 5))
        << "t=" << t;
  }
}

TEST_F(ЧебышёвPicardIteratorTest, Backwards) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = LinearProblem();
  Time const step = -1 * Second;  // Backwards!

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIterator<ODE> const integrator(ЧебышёвPicardIterationParams{
      .M = 64,
      .N = 64,
      .max_iterations = 64,
      .stopping_criterion = 2e-16,
  });

  auto const instance =
      integrator.NewInstance(problem.problem, append_state, step);
  auto const t_final = problem.t₀() + 2.5 * step;

  EXPECT_OK(instance->Solve(t_final));

  // Verify we have reached the desired final time.
  EXPECT_LE(states.back().s.value, t_final);  // ≤ because we are backwards.

  // Verify the results are close to the known solution.
  for (const auto& state : states) {
    auto t = state.s.value;
    auto y = std::get<0>(state.y).value;
    EXPECT_THAT(y, AlmostEquals(std::get<0>(problem.solution(t)), 0, 12))
        << "t=" << t;
  }
}

TEST_F(ЧебышёвPicardIteratorTest, Divergence) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = LinearProblem();
  Time const step = 60 * Second;  // Way too big; iteration won't converge.

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIterator<ODE> const integrator(ЧебышёвPicardIterationParams{
      .M = 64,
      .N = 64,
      .max_iterations = 64,
      .stopping_criterion = 0.1,  // Differences never even get this low!
  });

  auto const instance =
      integrator.NewInstance(problem.problem, append_state, step);
  auto const t_final = problem.t₀() + step;

  EXPECT_THAT(instance->Solve(t_final),
              StatusIs(absl::StatusCode::kFailedPrecondition));
}

TEST_F(ЧебышёвPicardIteratorTest, WriteToMessage) {
  ЧебышёвPicardIterator<ODE> const integrator(ЧебышёвPicardIterationParams{
      .M = 9,
      .N = 7,
      .max_iterations = 42,
      .stopping_criterion = 0.125,
  });

  serialization::FixedStepSizeIntegrator expected;
  expected.set_kind(serialization::FixedStepSizeIntegrator::CHEBYSHEV_PICARD);
  expected.mutable_chebyshev_picard_params()->set_m(9);
  expected.mutable_chebyshev_picard_params()->set_n(7);
  expected.mutable_chebyshev_picard_params()->set_max_iterations(42);
  expected.mutable_chebyshev_picard_params()->set_stopping_criterion(0.125);

  serialization::FixedStepSizeIntegrator actual;
  integrator.WriteToMessage(&actual);

  EXPECT_THAT(actual, EqualsProto(expected));
}

// Although in theory the iteration for y′ = y should converge for
// intervals
// of length < 40 for sufficiently high N (see [Bai10]), in
// practice the convergence degrades far earlier.
INSTANTIATE_TEST_SUITE_P(Linear, ЧебышёвPicardIteratorTest,
                         Values(
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-16,
                                 .step = 1 * Second,
                                 .tolerance = 9e-16,
                             },
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-16,
                                 .step = 2 * Second,
                                 .tolerance = 4.5e-15,
                             },
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-16,
                                 .step = 4 * Second,
                                 .tolerance = 8e-14,
                             },
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-12,
                                 .step = 8 * Second,
                                 .tolerance = 4e-10,
                             },
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = LinearProblem(),
                                 .N = 64,
                                 .stopping_criterion = 1e-7,
                                 .step = 16 * Second,
                                 // Absolute error is high due to the
                                 // exponential growth.
                                 .tolerance = 4e-3,
                             }));

// This problem appears in [BL69].
INSTANTIATE_TEST_SUITE_P(Tangent, ЧебышёвPicardIteratorTest,
                         Values(
                             // These are the parameters from [BL69].
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = TangentProblem(),
                                 .N = 16,
                                 .stopping_criterion = 0.5e-10,
                                 .step = 2 * Second,
                                 .tolerance = 8.7e-10,
                             },
                             // We can achieve better accuracy with higher N and
                             // a more stringent stopping criterion.
                             ЧебышёвPicardIteratorTestParam{
                                 .problem = TangentProblem(),
                                 .N = 32,
                                 .stopping_criterion = 1e-16,
                                 .step = 2 * Second,
                                 .tolerance = 4.5e-16,
                             }));

// From [Bai10]. Figures 9, 10, 12, 13 are slow and omitted for now.
INSTANTIATE_TEST_SUITE_P(
    PerturbedSinusoid, ЧебышёвPicardIteratorTest,
    Values(
        // [Bai10], figure 8.
        ЧебышёвPicardIteratorTestParam{
            .problem = PerturbedSinusoidProblem(/*ε=*/0.01,
                                                /*y₀=*/1),
            .N = 500,
            .stopping_criterion = 1e-16,
            .step = 64 * π * Second,
            .tolerance = 1e-9,
        },
        // [Bai10], figure 11.
        ЧебышёвPicardIteratorTestParam{
            .problem = PerturbedSinusoidProblem(/*ε=*/0.001,
                                                /*y₀=*/1),
            .N = 400,
            .stopping_criterion = 1e-16,
            .step = 64 * π * Second,
            .tolerance = 2.5e-11,
        }));

}  // namespace

}  // namespace integrators
}  // namespace principia
