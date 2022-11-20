#include "numerics/global_optimization.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
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

using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::Vector;
using quantities::Inverse;
using quantities::Length;
using quantities::Pow;
using quantities::si::Metre;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::AlmostEquals;
using testing_utilities::Branin;
using testing_utilities::Componentwise;
using testing_utilities::GoldsteinPrice;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::𝛁Branin;
using testing_utilities::𝛁GoldsteinPrice;
using testing_utilities::operator""_;
using ::testing::ElementsAre;
using ::testing::UnorderedElementsAre;

// The test functions in this file are from
// https://www.sfu.ca/~ssurjano/optimization.html.
class GlobalOptimizationTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
};

TEST_F(GlobalOptimizationTest, Branin) {
  using Optimizer = MultiLevelSingleLinkage<double, Displacement<World>>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto branin =
      [&function_invocations](Displacement<World> const& displacement) {
    ++function_invocations;
    auto const& coordinates = displacement.coordinates();
    // The extra |x₀| term ensures that we have a unique solution in three
    // dimensions.
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Pow<2>(x₀) + Branin(x₁, x₂);
  };

  auto grad_branin = [&gradient_invocations](
                         Displacement<World> const& displacement) {
    ++gradient_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 2 * x₀;
    auto const [g₁, g₂] = 𝛁Branin(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({0 * Metre, 2.5 * Metre, 7.5 * Metre}),
      .vertices = {
          Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 7.5 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 7.5 * Metre}),
      }};

  const auto tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, branin, grad_branin);
  auto const minima = optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                                 /*number_of_rounds=*/10,
                                                 tolerance);

  EXPECT_EQ(1434, function_invocations);
  EXPECT_EQ(598, gradient_invocations);

  // Note that the fourth minima is outside the |box| passed to the optimizer.
  EXPECT_THAT(
      minima,
      ElementsAre(
          Componentwise(
              AbsoluteErrorFrom(0 * Metre, IsNear(1.4e-7_(1) * Metre)),
              AbsoluteErrorFrom(9.42478 * Metre, IsNear(2.0e-6_(1) * Metre)),
              RelativeErrorFrom(2.475 * Metre, IsNear(8.0e-9_(1)))),
          Componentwise(
              AbsoluteErrorFrom(0 * Metre, IsNear(5.7e-7_(1) * Metre)),
              AbsoluteErrorFrom(π * Metre, IsNear(4.8e-9_(1) * Metre)),
              RelativeErrorFrom(2.275 * Metre, IsNear(9.1e-8_(1)))),
          Componentwise(
              AbsoluteErrorFrom(0 * Metre, IsNear(5.9e-8_(1) * Metre)),
              RelativeErrorFrom(-π * Metre, IsNear(3.5e-8_(1))),
              RelativeErrorFrom(12.275 * Metre, IsNear(6.2e-9_(1)))),
          Componentwise(
              AbsoluteErrorFrom(0 * Metre, IsNear(1.9e-8_(1) * Metre)),
              RelativeErrorFrom(5 * π * Metre, IsNear(7.3e-10_(1))),
              RelativeErrorFrom(12.875 * Metre, IsNear(1.0e-9_(1))))));
}

TEST_F(GlobalOptimizationTest, GoldsteinPrice) {
  using Optimizer = MultiLevelSingleLinkage<double, Displacement<World>>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto goldstein_price = [&function_invocations](
                             Displacement<World> const& displacement) {
    ++function_invocations;
    auto const& coordinates = displacement.coordinates();
    // The extra |x₀| term ensures that we have a unique solution in three
    // dimensions.
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Pow<2>(x₀) + GoldsteinPrice(x₁, x₂);
  };

  auto grad_goldstein_price = [&gradient_invocations](
                                  Displacement<World> const& displacement) {
    ++gradient_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 2 * x₀;
    auto const [g₁, g₂] = 𝛁GoldsteinPrice(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>(),
      .vertices = {
          Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
      }};

  const auto tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, goldstein_price, grad_goldstein_price);
  auto const minima = optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                                 /*number_of_rounds=*/10,
                                                 tolerance);

  EXPECT_EQ(2739, function_invocations);
  EXPECT_EQ(1812, gradient_invocations);
  EXPECT_THAT(minima,
              UnorderedElementsAre(
                  Componentwise(
                      AbsoluteErrorFrom(0 * Metre, IsNear(7.6e-7_(1) * Metre)),
                      AbsoluteErrorFrom(0 * Metre, IsNear(5.3e-8_(1) * Metre)),
                      RelativeErrorFrom(-1 * Metre, IsNear(3.8e-8_(1)))),
                  Componentwise(
                      AbsoluteErrorFrom(0 * Metre, IsNear(5.6e-8_(1) * Metre)),
                      RelativeErrorFrom(-0.6 * Metre, IsNear(6.8e-10_(1))),
                      RelativeErrorFrom(-0.4 * Metre, IsNear(1.1e-9_(1)))),
                  Componentwise(
                      AbsoluteErrorFrom(0 * Metre, IsNear(5.6e-8_(1) * Metre)),
                      RelativeErrorFrom(1.8 * Metre, IsNear(1.8e-10_(1))),
                      RelativeErrorFrom(0.2 * Metre, IsNear(7.0e-10_(1))))));
}

}  // namespace numerics
}  // namespace principia
