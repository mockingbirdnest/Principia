#if !_DEBUG

#include "integrators/чебышёв_picard_integrator.hpp"

#include <algorithm>
#include <concepts>
#include <vector>

#include "base/algebra.hpp"
#include "geometry/direct_sum.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/concepts.hpp"
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

using namespace principia::base::_algebra;
using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_ordinary_differential_equations;
using namespace principia::integrators::_чебышёв_picard_integrator;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

using FirstOrderODE =
    ExplicitFirstOrderOrdinaryDifferentialEquation<Instant, Position<World>>;
using SecondOrderODE =
    ExplicitSecondOrderOrdinaryDifferentialEquation<Position<World>>;

namespace {

// Returns a displacement of one metre along the x axis.
//
// Used to embed one-dimensional problems in the World.
Displacement<World> displacement() {
  return Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre});
}

template<quantity Scalar, typename Frame>
Scalar LInfinityNorm(Vector<Scalar, Frame> const& v) {
  return std::max(
      {Abs(v.coordinates().x), Abs(v.coordinates().y), Abs(v.coordinates().z)});
}

template<affine Scalar>
auto LInfinityNorm(DirectSum<Scalar> const& v) {
  return LInfinityNorm(get<0>(v));
}

template<typename Scalar>
auto LInfinityNorm(std::vector<Scalar> const& v) {
  auto norm = decltype(LInfinityNorm(std::declval<Scalar>()))();
  for (auto const& vᵢ : v) {
    norm = std::max(norm, LInfinityNorm(vᵢ));
  }
  return norm;
}

// An initial value problem with a known solution.
template<typename ODE_>
struct SolvedInitialValueProblem {
  using ODE = ODE_;
  using IndependentVariable = ODE::IndependentVariable;

  IndependentVariable t₀() const {
    return _чебышёв_picard_integrator::internal::
        ЧебышёвPicardIterationState<1, ODE>::GetIndependentVariable(
            problem.initial_state);
  }

  InitialValueProblem<ODE> problem;
  std::function<typename ODE::State(IndependentVariable const&)> solution;
};

// The first-order ODE y′ = y with y₀ = 1.
//
// The solution is y = eᵗ.
SolvedInitialValueProblem<FirstOrderODE> LinearProblem() {
  using ODE = FirstOrderODE;
  ODE linear_ode;
  linear_ode.compute_derivative =
      [](Instant const& t,
         ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = (y - World::origin) / Second;
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = linear_ode;
  problem.initial_state = {Instant(), {World::origin + displacement()}};

  return SolvedInitialValueProblem<ODE>{
      .problem = problem, .solution = [](Instant const& t) {
        return ODE::State(t,
                          std::exp((t - Instant()) / Second) * displacement() +
                              World::origin);
      }};
}

// The second-order ODE y″ = -y with y(0) = 0 and yʹ(0) = 1.
//
// The solution is y = sin(t).
SolvedInitialValueProblem<SecondOrderODE> SecondOrderLinearProblem() {
  using ODE = SecondOrderODE;
  ODE ode;
  ode.compute_acceleration =
      [](Instant const& t,
         ODE::DependentVariables const& positions,
         ODE::DependentVariableDerivatives const& velocities,
         ODE::DependentVariableDerivatives2& accelerations) {
        auto const& y = positions[0];
        auto& yʺ = accelerations[0];
        yʺ = -(y - World::origin) / (Second * Second);
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  problem.initial_state = {
      Instant(), {World::origin}, {displacement() / Second}};

  return SolvedInitialValueProblem<ODE>{
      .problem = problem, .solution = [](Instant const& t) -> ODE::State {
        Angle const θ = (t - Instant()) * Radian / Second;
        return ODE::State(t,
                          {Sin(θ) * displacement() + World::origin},
                          {Cos(θ) * displacement() / Second});
      }};
}

// The first-order ODE y′ = π/8 (1 + y²) with y(-1) = 0.
//
// The solution is y = tan(π/8 (t + 1)).
SolvedInitialValueProblem<FirstOrderODE> TangentProblem() {
  using ODE = FirstOrderODE;
  ODE ode;
  ode.compute_derivative =
      [](Instant const& t,
         ODE::DependentVariables const& dependent_variables,
         ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        double const y_unitless =
            InnerProduct(y - World::origin, displacement()) / (Metre * Metre);
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = π / 8 * (1 + y_unitless * y_unitless) / Second * displacement();
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  problem.initial_state = {Instant() - 1 * Second, {World::origin}};

  return SolvedInitialValueProblem<ODE>{
      .problem = problem, .solution = [](Instant const& t) {
        return ODE::State(t,
                          Tan(π / 8 * ((t - Instant()) / Second + 1) * Radian) *
                                  displacement() +
                              World::origin);
      }};
}

// The first-order ODE y′ = cos(t + εy).
//
// The solution ([Fuk97], equations 31 and 32) is
// y = -γt + 2/ε tan⁻¹((β + σ cos φ)/(1 + σ sin φ))
// where
// σ = α(sin φ + β cos φ)
// φ = ½ (1 - γε) t
// α = 2ε/(1 - ε + √(1 - ε²))
// β = 1/(1 + α) tan(εy₀/2)
// γ = ε/(1 + √(1 - ε²))
SolvedInitialValueProblem<FirstOrderODE> PerturbedSinusoidProblem(
    double const ε,
    double const y₀) {
  using ODE = FirstOrderODE;
  ODE ode;
  ode.compute_derivative =
      [ε](Instant const& t,
          ODE::DependentVariables const& dependent_variables,
          ODE::DependentVariableDerivatives& dependent_variable_derivatives) {
        auto const& [y] = dependent_variables;
        double const y_unitless =
            InnerProduct(y - World::origin, displacement()) / (Metre * Metre);
        auto& [yʹ] = dependent_variable_derivatives;
        yʹ = Cos(((t - Instant()) / Second + ε * y_unitless) * Radian) /
             Second * displacement();
        return absl::OkStatus();
      };

  InitialValueProblem<ODE> problem;
  problem.equation = ode;
  problem.initial_state = {Instant(), {y₀ * displacement() + World::origin}};

  double const one_plus_sqrt_one_minus_ε² = 1 + Sqrt(1 - ε * ε);
  double const α = 2 * ε / (one_plus_sqrt_one_minus_ε² - ε);
  double const β = 1 / (1 + α) * Tan(ε * y₀ / 2 * Radian);
  double const γ = ε / one_plus_sqrt_one_minus_ε²;

  return SolvedInitialValueProblem<ODE>{
      .problem = problem, .solution = [ε, α, β, γ](Instant const& t) {
        Angle const φ = 0.5 * (1 - γ * ε) * (t - Instant()) / Second * Radian;
        double const σ = α * (Sin(φ) + β * Cos(φ));

        return ODE::State(
            t,
            (-γ * (t - Instant()) / Second +
             2 / ε * ArcTan((β + σ * Cos(φ)) / (1 + σ * Sin(φ))) / Radian) *
                    displacement() +
                World::origin);
      }};
}

TEST(ЧебышёвPicardIntegratorTest, MultipleSteps) {
  // Set up the initial value problem.
  using ODE = FirstOrderODE;
  auto const problem = LinearProblem();
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
              [](auto const& Δ) { return LInfinityNorm(Δ) < 1e-16 * Metre; },
      });
  auto const t_final = problem.t₀() + 2.5 * step;

  EXPECT_OK(instance->Solve(t_final));

  // Verify we have reached the desired final time.
  EXPECT_GE(states.back().s.value, t_final);

  // Verify the results are close to the known solution.
  for (auto const& state : states) {
    auto t = state.s.value;
    auto y = get<0>(state.y).value;
    EXPECT_THAT(y, AlmostEquals(get<0>(problem.solution(t).y).value, 0, 5))
        << "t=" << t;
  }
}

TEST(ЧебышёвPicardIntegratorTest, Backwards) {
  // Set up the initial value problem.
  using ODE = FirstOrderODE;
  auto const problem = LinearProblem();
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
              [](auto const& Δ) { return LInfinityNorm(Δ) < 2e-16 * Metre; },
      });
  auto const t_final = problem.t₀() + 2.5 * step;

  EXPECT_OK(instance->Solve(t_final));

  // Verify we have reached the desired final time.
  EXPECT_LE(states.back().s.value, t_final);  // ≤ because we are backwards.

  // Verify the results are close to the known solution.
  for (auto const& state : states) {
    auto t = state.s.value;
    auto y = get<0>(state.y).value;
    EXPECT_THAT(y, AlmostEquals(get<0>(problem.solution(t).y).value, 0, 12))
        << "t=" << t;
  }
}

TEST(ЧебышёвPicardIntegratorTest, Divergence) {
  // Set up the initial value problem.
  using ODE = FirstOrderODE;
  auto const problem = LinearProblem();
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
                                       return LInfinityNorm(Δ) <= 0.1 * Metre;
                                     },
                             });
  auto const t_final = problem.t₀() + step;

  EXPECT_THAT(instance->Solve(t_final),
              StatusIs(absl::StatusCode::kFailedPrecondition));
}

template<typename T>
concept ЧебышёвPicardIntegratorTestParameters = requires {
  // The ODE type, which should be FirstOrderODE or SecondOrderODE.
  typename T::ODE;

  // The ODE (with solution) under test.
  {
    T::problem()
  } -> std::convertible_to<SolvedInitialValueProblem<typename T::ODE>>;

  // The ЧебышёвPicard method.
  typename T::Method;

  // The step size used for this test.
  {
    T::step
  } -> std::convertible_to<typename T::ODE::IndependentVariableDifference>;

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
  using ODE = typename TypeParam::ODE;
  auto const problem = TypeParam::problem();
  Time const step = TypeParam::step;

  std::vector<typename ODE::State> states;
  auto const append_state = [&states](typename ODE::State const& state) {
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
                if constexpr (ODE::order == 1) {
                  return LInfinityNorm(Δ) < TypeParam::stopping_criterion;
                } else {
                  return LInfinityNorm(Δ.position_error) <
                             TypeParam::stopping_criterion &&
                         LInfinityNorm(Δ.velocity_error) <
                             TypeParam::stopping_criterion / Second;
                }
              },
      });
  EXPECT_OK(instance->Solve(problem.t₀() + step));

  // Verify the results are close to the known solution.
  for (auto const& state : states) {
    if constexpr (ODE::order == 1) {
      auto const t = state.s.value;
      auto const y = get<0>(state.y).value;
      EXPECT_THAT(y,
                  AbsoluteErrorFrom(get<0>(problem.solution(t).y).value,
                                    Lt(TypeParam::tolerance)))
          << "t=" << t;
    } else {  // ODE::order == 2
      auto const t = state.time.value;
      auto const actual_q = state.positions[0].value;
      auto const actual_v = state.velocities[0].value;

      auto const solution = problem.solution(t);
      auto const expected_q = solution.positions[0].value;
      auto const expected_v = solution.velocities[0].value;
      EXPECT_THAT(actual_q,
                  AbsoluteErrorFrom(expected_q, Lt(TypeParam::tolerance)))
          << "t=" << t;
      EXPECT_THAT(
          actual_v,
          AbsoluteErrorFrom(expected_v, Lt(TypeParam::tolerance / Second)))
          << "t=" << t;
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(ЧебышёвPicardIntegratorParameterizedTest,
                            Convergence);

// Although in theory the iteration for y′ = y should converge for
// intervals
// of length < 40 for sufficiently high N (see [BJ12]), in
// practice the convergence degrades far earlier.
struct Linear1Second : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 1 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 9e-16 * Metre;
};
struct Linear2Seconds : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 2 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 4.5e-15 * Metre;
};
struct Linear4Seconds : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 4 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 8e-14 * Metre;
};
struct Linear8Seconds : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return LinearProblem();
  }

  using Method = ЧебышёвPicard<64>;
  static constexpr Time step = 8 * Second;
  static constexpr Length stopping_criterion = 1e-12 * Metre;
  static constexpr Length tolerance = 4e-10 * Metre;
};
struct Linear16Seconds : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
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
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
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

struct SecondOrderLinear1Second : not_constructible {
  using ODE = SecondOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return SecondOrderLinearProblem();
  }

  using Method = ЧебышёвPicard<32>;
  static constexpr Time step = 1 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 3e-16 * Metre;
};
struct SecondOrderLinear2Seconds : not_constructible {
  using ODE = SecondOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return SecondOrderLinearProblem();
  }

  using Method = ЧебышёвPicard<32>;
  static constexpr Time step = 2 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 5e-16 * Metre;
};
struct SecondOrderLinear4Seconds : not_constructible {
  using ODE = SecondOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return SecondOrderLinearProblem();
  }

  using Method = ЧебышёвPicard<32>;
  static constexpr Time step = 4 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 2e-15 * Metre;
};
struct SecondOrderLinear8Seconds : not_constructible {
  using ODE = SecondOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return SecondOrderLinearProblem();
  }

  using Method = ЧебышёвPicard<32>;
  static constexpr Time step = 8 * Second;
  static constexpr Length stopping_criterion = 1e-12 * Metre;
  static constexpr Length tolerance = 3e-14 * Metre;
};
struct SecondOrderLinearMGreaterThanN : not_constructible {
  using ODE = SecondOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return SecondOrderLinearProblem();
  }

  using Method = ЧебышёвPicard<64, 32>;
  static constexpr Time step = 1 * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 3e-16 * Metre;
};

using SecondOrderLinear = Types<SecondOrderLinear1Second,
                                SecondOrderLinear2Seconds,
                                SecondOrderLinear4Seconds,
                                SecondOrderLinear8Seconds,
                                SecondOrderLinearMGreaterThanN>;
INSTANTIATE_TYPED_TEST_SUITE_P(SecondOrderLinear,
                               ЧебышёвPicardIntegratorParameterizedTest,
                               SecondOrderLinear);

// This problem appears in [BL69].
struct TangentA : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return TangentProblem();
  }

  // These are the parameters from [BL69].
  using Method = ЧебышёвPicard<16>;
  static constexpr Time step = 2 * Second;
  static constexpr Length stopping_criterion = 0.5e-10 * Metre;
  static constexpr Length tolerance = 8.7e-10 * Metre;
};
struct TangentB : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
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
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
    return PerturbedSinusoidProblem(/*ε=*/0.01,
                                    /*y₀=*/1);
  }

  using Method = ЧебышёвPicard<500>;
  static constexpr Time step = 64 * π * Second;
  static constexpr Length stopping_criterion = 1e-16 * Metre;
  static constexpr Length tolerance = 1e-9 * Metre;
};
struct PerturbedSinusoidFigure6 : not_constructible {
  using ODE = FirstOrderODE;
  static SolvedInitialValueProblem<ODE> problem() {
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
