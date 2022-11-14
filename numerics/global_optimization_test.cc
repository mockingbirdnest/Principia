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
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_;
using ::testing::Lt;
using ::testing::UnorderedElementsAre;

// The test functions in this file are from
// https://www.sfu.ca/~ssurjano/optimization.html.
class GlobalOptimizationTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
};

TEST_F(GlobalOptimizationTest, Smoke) {
  using Foo = MultiLevelSingleLinkage<Length, Position<World>>;
}

TEST_F(GlobalOptimizationTest, GoldsteinPrice) {
  using Optimizer = MultiLevelSingleLinkage<double, Displacement<World>>;
  int function_invocations = 0;
  int gradient_invocations = 0;

  auto goldstein_price = [&function_invocations](
                             Displacement<World> const& displacement) {
    ++function_invocations;
    auto const& coordinates = displacement.coordinates();
    // The extra |x0| term ensures that we have a unique solution in three
    // dimensions.
    double const x0 = coordinates[0] / Metre;
    double const x1 = coordinates[1] / Metre;
    double const x2 = coordinates[2] / Metre;
    return Pow<2>(x0) +
           (1 + Pow<2>(x1 + x2 + 1) * (19 - 14 * x1 + 3 * Pow<2>(x1) - 14 * x2 +
                                       6 * x1 * x2 + 3 * Pow<2>(x2))) *
               (30 + Pow<2>(2 * x1 - 3 * x2) *
                         (18 - 32 * x1 + 12 * Pow<2>(x1) + 48 * x2 -
                          36 * x1 * x2 + 27 * Pow<2>(x2)));
  };

  auto grad_goldstein_price = [&gradient_invocations](
                                  Displacement<World> const& displacement) {
    ++gradient_invocations;
    auto const& coordinates = displacement.coordinates();
    double const x0 = coordinates[0] / Metre;
    double const x1 = coordinates[1] / Metre;
    double const x2 = coordinates[2] / Metre;
    double const g0 = 2 * x0;
    double const g1 =
        24 * (-1 + 2 * x1 - 3 * x2) * (2 * x1 - 3 * x2) *
            (2 * x1 - 3 * (1 + x2)) *
            (1 +
             Pow<2>(1 + x1 + x2) * (19 + 3 * Pow<2>(x1) + x2 * (-14 + 3 * x2) +
                                    2 * x1 * (-7 + 3 * x2))) +
        12 * (-2 + x1 + x2) * (-1 + x1 + x2) * (1 + x1 + x2) *
            (30 + Pow<2>(2 * x1 - 3 * x2) *
                      (18 + 12 * Pow<2>(x1) - 4 * x1 * (8 + 9 * x2) +
                       3 * x2 * (16 + 9 * x2)));
    double const g2 =
        -36 * (-1 + 2 * x1 - 3 * x2) * (2 * x1 - 3 * x2) *
            (2 * x1 - 3 * (1 + x2)) *
            (1 +
             Pow<2>(1 + x1 + x2) * (19 + 3 * Pow<2>(x1) + x2 * (-14 + 3 * x2) +
                                    2 * x1 * (-7 + 3 * x2))) +
        12 * (-2 + x1 + x2) * (-1 + x1 + x2) * (1 + x1 + x2) *
            (30 + Pow<2>(2 * x1 - 3 * x2) *
                      (18 + 12 * Pow<2>(x1) - 4 * x1 * (8 + 9 * x2) +
                       3 * x2 * (16 + 9 * x2)));
    return Vector<Inverse<Length>, World>({g0 / Metre, g1 / Metre, g2 / Metre});
  };

  // Correctness checks for the function and its gradient.
  {
    auto const test_point =
        Displacement<World>({0 * Metre, 0.5 * Metre, -0.3 * Metre});
    EXPECT_THAT(goldstein_price(test_point), IsNear(596.161_(1)));
    EXPECT_THAT(grad_goldstein_price(test_point),
                Componentwise(AlmostEquals(0 / Metre, 0),
                              IsNear(-601.51_(1) / Metre),
                              IsNear(2163.65_(1) / Metre)));
  }

  Optimizer::Box const box = {
      .centre = Displacement<World>(),
      .vertices = {
          Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
      }};

  const auto tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, goldstein_price, grad_goldstein_price);
  auto const minima = optimizer.FindGlobalMinima(10, 10, tolerance);

  EXPECT_EQ(32278, function_invocations);
  EXPECT_EQ(17293, gradient_invocations);
  EXPECT_THAT(
      minima,
      UnorderedElementsAre(
          Componentwise(AbsoluteErrorFrom(0 * Metre, Lt(9.9e-8 * Metre)),
                        AbsoluteErrorFrom(0 * Metre, Lt(1.5e-7 * Metre)),
                        RelativeErrorFrom(-1 * Metre, Lt(8.0e-8))),
          Componentwise(AbsoluteErrorFrom(0 * Metre, Lt(5.6e-8 * Metre)),
                        RelativeErrorFrom(-0.6 * Metre, Lt(6.8e-10)),
                        RelativeErrorFrom(-0.4 * Metre, Lt(1.1e-9))),
          Componentwise(AbsoluteErrorFrom(0 * Metre, Lt(9.4e-8 * Metre)),
                        RelativeErrorFrom(1.2 * Metre, Lt(2.8e-9)),
                        RelativeErrorFrom(0.8 * Metre, Lt(3.0e-9))),
          Componentwise(AbsoluteErrorFrom(0 * Metre, Lt(4.9e-7 * Metre)),
                        RelativeErrorFrom(1.8 * Metre, Lt(1.4e-7)),
                        RelativeErrorFrom(0.2 * Metre, Lt(5.0e-7)))));
}

}  // namespace numerics
}  // namespace principia
