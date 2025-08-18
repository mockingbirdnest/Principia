#include "numerics/global_optimization.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "quantities/arithmetic.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/optimization_test_functions.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace numerics {

using ::testing::AnyOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::_;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_space;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_global_optimization;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;
using namespace principia::testing_utilities::_optimization_test_functions;
using namespace principia::testing_utilities::_vanishes_before;

// The test functions in this file are from
// https://www.sfu.ca/~ssurjano/optimization.html.
class GlobalOptimizationTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag>;
};

TEST_F(GlobalOptimizationTest, Branin) {
  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/2>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto branin =
      [&function_invocations](Displacement<World> const& displacement) {
    ++function_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Branin(x₁, x₂);
  };

  auto grad_branin = [&gradient_invocations](
                         Displacement<World> const& displacement) {
    ++gradient_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 0;
    auto const [g₁, g₂] = 𝛁Branin(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({0 * Metre, 2.5 * Metre, 7.5 * Metre}),
      .vertices = {
          Displacement<World>({0 * Metre, 7.5 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 7.5 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, branin, grad_branin);
  {
    auto const minima = optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                                   /*number_of_rounds=*/10,
                                                   tolerance);

    EXPECT_THAT(function_invocations,
                AnyOf(Eq(1264),    // MSVC.
                      Eq(1243)));  // Clang.
    EXPECT_EQ(316, gradient_invocations);

    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                _,
                AbsoluteErrorFrom(9.42478 * Metre, IsNear(1.8e-6_(1) * Metre)),
                AbsoluteErrorFrom(2.475 * Metre, IsNear(2.0e-7_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(-π * Metre, IsNear(2.8e-7_(1) * Metre)),
                AbsoluteErrorFrom(12.275 * Metre, IsNear(4.9e-8_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(π * Metre, IsNear(3.3e-8_(1) * Metre)),
                AbsoluteErrorFrom(2.275 * Metre, IsNear(5.7e-7_(1) * Metre)))));
  }
  function_invocations = 0;
  gradient_invocations = 0;
  {
    auto const minima =
        optimizer.FindGlobalMinima(/*points_per_round=*/100,
                                   /*number_of_rounds=*/std::nullopt,
                                   tolerance);

    EXPECT_THAT(function_invocations,
                AnyOf(Eq(1019),   // MSVC.
                      Eq(997)));  // Clang.
    EXPECT_EQ(136, gradient_invocations);

    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                _,
                AbsoluteErrorFrom(π * Metre, IsNear(6.9e-9_(1) * Metre)),
                AbsoluteErrorFrom(2.275 * Metre, IsNear(6.9e-9_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(-π * Metre, IsNear(4.8e-8_(1) * Metre)),
                AbsoluteErrorFrom(12.275 * Metre, IsNear(1.5e-8_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(9.42478 * Metre, IsNear(1.6e-6_(1) * Metre)),
                AbsoluteErrorFrom(2.475 * Metre, IsNear(2.9e-7_(1) * Metre)))));
  }
}

TEST_F(GlobalOptimizationTest, GoldsteinPrice) {
  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/2>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto goldstein_price = [&function_invocations](
                             Displacement<World> const& displacement) {
    ++function_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return GoldsteinPrice(x₁, x₂);
  };

  auto grad_goldstein_price = [&gradient_invocations](
                                  Displacement<World> const& displacement) {
    ++gradient_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 0;
    auto const [g₁, g₂] = 𝛁GoldsteinPrice(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>(),
      .vertices = {
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, goldstein_price, grad_goldstein_price);

  {
    auto const minima = optimizer.FindGlobalMinima(/*points_per_round=*/20,
                                                   /*number_of_rounds=*/10,
                                                   tolerance);

    EXPECT_THAT(function_invocations,
                AnyOf(Eq(2129),    // MSVC.
                      Eq(2078)));  // Clang.
    EXPECT_EQ(278, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                _,
                AbsoluteErrorFrom(-0.6 * Metre, IsNear(9.6e-8_(1) * Metre)),
                AbsoluteErrorFrom(-0.4 * Metre, IsNear(5.8e-9_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(0 * Metre, IsNear(6.2e-8_(1) * Metre)),
                AbsoluteErrorFrom(-1 * Metre, IsNear(1.1e-7_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(1.8 * Metre, IsNear(5.6e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.2 * Metre, IsNear(4.8e-7_(1) * Metre)))));
  }
  function_invocations = 0;
  gradient_invocations = 0;
  {
    auto const minima =
        optimizer.FindGlobalMinima(/*points_per_round=*/500,
                                   /*number_of_rounds=*/std::nullopt,
                                   tolerance);

    EXPECT_THAT(function_invocations,
                AnyOf(Eq(5730),    // MSVC.
                      Eq(5674)));  // Clang.
    EXPECT_EQ(178, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                _,
                AbsoluteErrorFrom(0 * Metre, IsNear(4.4e-8_(1) * Metre)),
                AbsoluteErrorFrom(-1 * Metre, IsNear(5.9e-10_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(1.8 * Metre, IsNear(2.8e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.2 * Metre, IsNear(5.8e-7_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(-0.6 * Metre, IsNear(3.3e-8_(1) * Metre)),
                AbsoluteErrorFrom(-0.4 * Metre, IsNear(6.1e-8_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(1.2 * Metre, IsNear(5.3e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.8 * Metre, IsNear(4.3e-7_(1) * Metre)))));
  }
}

TEST_F(GlobalOptimizationTest, Hartmann3) {
  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/3>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto hartmann3 =
      [&function_invocations](Displacement<World> const& displacement) {
    ++function_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Hartmann3(x₀, x₁, x₂);
  };

  auto grad_hartmann3 = [&gradient_invocations](
                            Displacement<World> const& displacement) {
    ++gradient_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    auto const [g₀, g₁, g₂] = 𝛁Hartmann3(x₀, x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({0.5 * Metre, 0.5 * Metre, 0.5 * Metre}),
      .vertices = {
          Displacement<World>({0.5 * Metre, 0 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0.5 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 0.5 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, hartmann3, grad_hartmann3);

  {
    auto const minima = optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                                   /*number_of_rounds=*/10,
                                                   tolerance);

    EXPECT_THAT(function_invocations,
                AnyOf(Eq(1319),    // MSVC.
                      Eq(1357)));  // Clang.
    EXPECT_EQ(429, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                AbsoluteErrorFrom(0.114614 * Metre, IsNear(2.0e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.555649 * Metre, IsNear(9.4e-8_(1) * Metre)),
                AbsoluteErrorFrom(0.852547 * Metre,
                                  IsNear(6.0e-8_(1) * Metre))),
            Componentwise(
                AbsoluteErrorFrom(0.109338 * Metre, IsNear(1.8e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.860524 * Metre, IsNear(1.8e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.564123 * Metre,
                                  IsNear(1.6e-7_(1) * Metre))),
            Componentwise(
                AbsoluteErrorFrom(0.368723 * Metre, IsNear(1.7e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.117561 * Metre, IsNear(4.8e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.267573 * Metre,
                                  IsNear(7.5e-7_(1) * Metre)))));
  }
  function_invocations = 0;
  gradient_invocations = 0;
  {
    auto const minima =
        optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                   /*number_of_rounds=*/std::nullopt,
                                   tolerance);

    EXPECT_EQ(199, function_invocations);
    EXPECT_EQ(126, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                AbsoluteErrorFrom(0.368723 * Metre, IsNear(2.1e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.117561 * Metre, IsNear(6.5e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.267573 * Metre,
                                  IsNear(7.3e-7_(1) * Metre))),
            Componentwise(
                AbsoluteErrorFrom(0.114614 * Metre, IsNear(1.6e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.555649 * Metre, IsNear(2.7e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.852547 * Metre,
                                  IsNear(2.1e-8_(1) * Metre)))));
  }
}

// A function that looks like the opposite of the gravitational potential.
TEST_F(GlobalOptimizationTest, Potential) {
  using Optimizer = MultiLevelSingleLinkage<Inverse<Length>,
                                            Displacement<World>,
                                            /*dimensions=*/3>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto potential =
      [&function_invocations](Displacement<World> const& displacement) {
    ++function_invocations;
    return 1 / displacement.Norm();
  };

  auto grad_potential =
      [&gradient_invocations](Displacement<World> const& displacement) {
    ++gradient_invocations;
    return -displacement / Pow<3>(displacement.Norm());
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({1 * Metre, 1 * Metre, 1 * Metre}),
      .vertices = {
          Displacement<World>({10 * Metre, 0 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 10 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 10 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, potential, grad_potential);

  {
    auto const minima = optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                                   /*number_of_rounds=*/10,
                                                   tolerance);

    EXPECT_THAT(function_invocations,
                AnyOf(Eq(1452),    // MSVC.
                      Eq(1448)));  // Clang.
    EXPECT_EQ(503, gradient_invocations);
    EXPECT_THAT(minima, IsEmpty());
  }
  function_invocations = 0;
  gradient_invocations = 0;
  {
    auto const minima =
        optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                   /*number_of_rounds=*/std::nullopt,
                                   tolerance);

    EXPECT_EQ(127, function_invocations);
    EXPECT_EQ(91, gradient_invocations);
    EXPECT_THAT(minima, IsEmpty());
  }
}

}  // namespace numerics
}  // namespace principia
