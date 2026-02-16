#if !_DEBUG

#include "integrators/чебышёв_picard_integrator.hpp"

#include <concepts>
#include <tuple>
#include <vector>

#include "geometry/instant.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace integrators {

using ::testing::Lt;
using ::testing::Test;
using ::testing::Types;
using ::testing::Values;

using namespace principia::geometry::_instant;
using namespace principia::geometry::_point;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_чебышёв_picard_integrator;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

using ODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Point<Length>>;

namespace {

// An initial value problem with a known solution.
struct SolvedInitialValueProblem {
  ODE::IndependentVariable t₀() const {
    return problem.initial_state.s.value;
  }

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
      [](Instant const& t,
         ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = (y - Point<Length>()) / Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = linear_ode;
  problem.initial_state = {Instant(), {Point<Length>() + 1 * Metre}};

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [](Instant const& t) -> std::tuple<Point<Length>> {
        return std::exp((t - Instant()) / Second) * Metre + Point<Length>();
      }};
}

// The first order ODE y′ = π/8 (1 + y²) with y(-1) = 0.
//
// The solution is y = tan(π/8 (t + 1)).
SolvedInitialValueProblem TangentProblem() {
  ODE ode;
  ode.compute_derivative =
      [](Instant const& t,
         ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        double const y_unitless = (y - Point<Length>()) / Metre;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = π / 8 * (1 + y_unitless * y_unitless) * Metre / Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  problem.initial_state = {Instant() - 1 * Second, {Point<Length>()}};

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [](Instant const& t) -> std::tuple<Point<Length>> {
        return Tan(π / 8 * ((t - Instant()) / Second + 1) * Radian) * Metre +
               Point<Length>();
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
SolvedInitialValueProblem PerturbedSinusoidProblem(double const ε,
                                                   double const y₀) {
  ODE ode;
  ode.compute_derivative =
      [ε](Instant const& t,
          ODE::DependentVariables const& dependent_variables,
          ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        double const y_unitless = (y - Point<Length>()) / Metre;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = Cos(((t - Instant()) / Second + ε * y_unitless) * Radian) * Metre /
             Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  problem.initial_state = {Instant(), {y₀ * Metre + Point<Length>()}};

  double const one_plus_sqrt_one_minus_ε² = 1 + Sqrt(1 - ε * ε);
  double const α = 2 * ε / (one_plus_sqrt_one_minus_ε² - ε);
  double const β = 1 / (1 + α) * Tan(ε * y₀ / 2 * Radian);
  double const γ = ε / one_plus_sqrt_one_minus_ε²;

  return SolvedInitialValueProblem{
      .problem = problem,
      .solution = [ε, α, β, γ](Instant const& t) -> std::tuple<Point<Length>> {
        Angle const φ = 0.5 * (1 - γ * ε) * (t - Instant()) / Second * Radian;
        double const σ = α * (Sin(φ) + β * Cos(φ));

        return (-γ * (t - Instant()) / Second +
                2 / ε * ArcTan((β + σ * Cos(φ)) / (1 + σ * Sin(φ))) / Radian) *
                   Metre +
               Point<Length>();
      }};
}

TEST(ЧебышёвPicardIntegratorTest, MultipleSteps) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = LinearProblem();
  Time const step = 1 * Second;

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIntegrator<ЧебышёвPicard<64, 64>, ODE> const integrator;

  auto const instance = integrator.NewInstance(
      problem.problem,
      append_state,
      step,
      ЧебышёвPicardIterationParams<ODE>{
          .max_iterations = 64,
          .stopping_criterion =
              [](auto const& Δ) { return std::get<0>(Δ) < 1e-16 * Metre; },
      });
  auto const t_final = problem.t₀() + 2.5 * step;

  EXPECT_OK(instance->Solve(t_final));

  // Verify we have reached the desired final time.
  EXPECT_GE(states.back().s.value, t_final);

  // Verify the results are close to the known solution.
  for (auto const& state : states) {
    auto t = state.s.value;
    auto y = std::get<0>(state.y).value;
    EXPECT_THAT(y, AlmostEquals(std::get<0>(problem.solution(t)), 0, 5))
        << "t=" << t;
  }
}

TEST(ЧебышёвPicardIntegratorTest, Backwards) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = LinearProblem();
  Time const step = -1 * Second;  // Backwards!

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIntegrator<ЧебышёвPicard<64, 64>, ODE> const integrator;

  auto const instance = integrator.NewInstance(
      problem.problem,
      append_state,
      step,
      ЧебышёвPicardIterationParams<ODE>{
          .max_iterations = 64,
          .stopping_criterion =
              [](auto const& Δ) { return std::get<0>(Δ) < 2e-16 * Metre; },
      });
  auto const t_final = problem.t₀() + 2.5 * step;

  EXPECT_OK(instance->Solve(t_final));

  // Verify we have reached the desired final time.
  EXPECT_LE(states.back().s.value, t_final);  // ≤ because we are backwards.

  // Verify the results are close to the known solution.
  for (auto const& state : states) {
    auto t = state.s.value;
    auto y = std::get<0>(state.y).value;
    EXPECT_THAT(y, AlmostEquals(std::get<0>(problem.solution(t)), 0, 12))
        << "t=" << t;
  }
}

TEST(ЧебышёвPicardIntegratorTest, Divergence) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const& problem = LinearProblem();
  Time const step = 60 * Second;  // Way too big; iteration won't converge.

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIntegrator<ЧебышёвPicard<64, 64>, ODE> const integrator;

  auto const instance =
      integrator.NewInstance(problem.problem,
                             append_state,
                             step,
                             ЧебышёвPicardIterationParams<ODE>{
                                 .max_iterations = 64,
                                 .stopping_criterion =
                                     [](auto const& Δ) {
                                       // Differences never even get this low!
                                       return std::get<0>(Δ) <= 0.1 * Metre;
                                     },
                             });
  auto const t_final = problem.t₀() + step;

  EXPECT_THAT(instance->Solve(t_final),
              StatusIs(absl::StatusCode::kFailedPrecondition));
}

template<typename T>
concept ЧебышёвPicardIntegratorTestParameters = requires {
  // The ODE (with solution) under test.
  { T::problem() } -> std::convertible_to<SolvedInitialValueProblem>;

  // The ЧебышёвPicard method.
  typename T::Method;

  // The step size used for this test.
  { T::step } -> std::convertible_to<ODE::IndependentVariableDifference>;

  // The stopping criterion used for this test.
  { T::stopping_criterion } -> std::convertible_to<Length>;

  // The distance the integrated solution is allowed to vary from the analytic
  // solution.
  { T::tolerance } -> std::convertible_to<Length>;
};

template<ЧебышёвPicardIntegratorTestParameters T>
class ЧебышёвPicardIntegratorParameterizedTest : public Test {};

TYPED_TEST_SUITE_P(ЧебышёвPicardIntegratorParameterizedTest);

TYPED_TEST_P(ЧебышёвPicardIntegratorParameterizedTest, Convergence) {
  // Set up the initial value problem.
  SolvedInitialValueProblem const problem = TypeParam::problem();
  Time const step = TypeParam::step;

  std::vector<ODE::State> states;
  auto const append_state = [&states](ODE::State const& state) {
    states.push_back(state);
  };

  // Build the integrator and solve the problem.
  ЧебышёвPicardIntegrator<typename TypeParam::Method, ODE> const integrator;

  auto const instance = integrator.NewInstance(
      problem.problem,
      append_state,
      step,
      ЧебышёвPicardIterationParams<ODE>{
          .max_iterations = 64,
          .stopping_criterion =
              [](auto const& Δ) {
                return std::get<0>(Δ) < TypeParam::stopping_criterion;
              },
      });
  EXPECT_OK(instance->Solve(problem.t₀() + step));

  // Verify the results are close to the known solution.
  for (auto const& state : states) {
    auto t = state.s.value;
    auto y = std::get<0>(state.y).value;
    EXPECT_THAT(y,
                AbsoluteErrorFrom(std::get<0>(problem.solution(t)),
                                  Lt(TypeParam::tolerance)))
        << "t=" << t;
  }
}

REGISTER_TYPED_TEST_SUITE_P(ЧебышёвPicardIntegratorParameterizedTest,
                            Convergence);

// Although in theory the iteration for y′ = y should converge for
// intervals
// of length < 40 for sufficiently high N (see [BJ12]), in
// practice the convergence degrades far earlier.
struct Linear1Second : not_constructible {
  static SolvedInitialValueProblem problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 1 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 9e-16 * Metre;
};
struct Linear2Seconds : not_constructible {
  static SolvedInitialValueProblem problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 2 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 4.5e-15 * Metre;
};
struct Linear4Seconds : not_constructible {
  static SolvedInitialValueProblem problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 4 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 8e-14 * Metre;
};
struct Linear8Seconds : not_constructible {
  static SolvedInitialValueProblem problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 8 * Second;
  static constexpr Length stopping_criterion = 1e-12 * Metre;
  static constexpr Length tolerance = 4e-10 * Metre;
};
struct Linear16Seconds : not_constructible {
  static SolvedInitialValueProblem problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 16 * Second;
  static constexpr Length stopping_criterion = 1e-7 * Metre;
  // Absolute error is high due to the
  // exponential growth.
  static constexpr Length tolerance = 4e-3 * Metre;
};
struct LinearMGreaterThanN : not_constructible {
  static SolvedInitialValueProblem problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<128, 64>;
  static constexpr Time step = 1 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 2e-15 * Metre;
};

using Linear = Types<Linear1Second,
                     Linear2Seconds,
                     Linear4Seconds,
                     Linear8Seconds,
                     Linear16Seconds,
                     LinearMGreaterThanN>;
INSTANTIATE_TYPED_TEST_SUITE_P(Linear,
                               ЧебышёвPicardIntegratorParameterizedTest,
                               Linear);

// This problem appears in [BL69].
struct TangentA : not_constructible {
  static SolvedInitialValueProblem problem() {
    return TangentProblem();
  }

  // These are the parameters from [BL69].
  using Method = ЧебышёвPicard<16>;
  static constexpr Time step = 2 * Second;
  static constexpr Length stopping_criterion = 0.5e-10 * Metre;
  static constexpr Length tolerance = 8.7e-10 * Metre;
};
struct TangentB : not_constructible {
  static SolvedInitialValueProblem problem() {
    return TangentProblem();
  }

  // We can achieve better accuracy with higher N and
  // a more stringent stopping criterion.
  using Method = ЧебышёвPicard<32>;
  static constexpr Time step = 2 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 4.5e-16 * Metre;
};

using Tangent = Types<TangentA, TangentB>;
INSTANTIATE_TYPED_TEST_SUITE_P(Tangent,
                               ЧебышёвPicardIntegratorParameterizedTest,
                               Tangent);

// From [BJ10]. Figures 4, 5, 7, 8 are slow and omitted for now.
struct PerturbedSinusoidFigure3 : not_constructible {
  static SolvedInitialValueProblem problem() {
    return PerturbedSinusoidProblem(/*ε=*/0.01,
                                    /*y₀=*/1);
  }

  using Method = ЧебышёвPicard<500>;
  static constexpr Time step = 64 * π * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 1e-9 * Metre;
};
struct PerturbedSinusoidFigure6 : not_constructible {
  static SolvedInitialValueProblem problem() {
    return PerturbedSinusoidProblem(/*ε=*/0.001,
                                    /*y₀=*/1);
  }

  using Method = ЧебышёвPicard<400>;
  static constexpr Time step = 64 * π * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 2.5e-11 * Metre;
};

using PerturbedSinusoid =
    Types<PerturbedSinusoidFigure3, PerturbedSinusoidFigure6>;
INSTANTIATE_TYPED_TEST_SUITE_P(PerturbedSinusoid,
                               ЧебышёвPicardIntegratorParameterizedTest,
                               PerturbedSinusoid);

}  // namespace

}  // namespace integrators
}  // namespace principia

#endif
