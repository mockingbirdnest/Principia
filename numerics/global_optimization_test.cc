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
using testing_utilities::Componentwise;
using testing_utilities::GradGoldsteinPrice;
using testing_utilities::GoldsteinPrice;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_;
using ::testing::UnorderedElementsAre;

// The test functions in this file are from
// https://www.sfu.ca/~ssurjano/optimization.html.
class GlobalOptimizationTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
};

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
    auto const [g₁, g₂] = GradGoldsteinPrice(x₁, x₂);
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
