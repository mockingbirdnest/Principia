#include "numerics/global_optimization.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
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
using testing_utilities::Hartmann3;
using testing_utilities::IsNear;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::𝛁Branin;
using testing_utilities::𝛁GoldsteinPrice;
using testing_utilities::𝛁Hartmann3;
using testing_utilities::operator""_;
using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::_;

// The test functions in this file are from
// https://www.sfu.ca/~ssurjano/optimization.html.
class GlobalOptimizationTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
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

    EXPECT_EQ(1264, function_invocations);
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

    EXPECT_EQ(1019, function_invocations);
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

    EXPECT_EQ(2131, function_invocations);
    EXPECT_EQ(280, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                _,
                AbsoluteErrorFrom(-0.6 * Metre, IsNear(9.6e-8_(1) * Metre)),
                AbsoluteErrorFrom(-0.4 * Metre, IsNear(5.8e-9_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(0 * Metre, IsNear(7.1e-8_(1) * Metre)),
                AbsoluteErrorFrom(-1 * Metre, IsNear(2.4e-7_(1) * Metre))),
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

    EXPECT_EQ(5738, function_invocations);
    EXPECT_EQ(187, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                _,
                AbsoluteErrorFrom(0 * Metre, IsNear(4.4e-8_(1) * Metre)),
                AbsoluteErrorFrom(-1 * Metre, IsNear(5.9e-10_(1) * Metre))),
            Componentwise(
                _,
                AbsoluteErrorFrom(1.8 * Metre, IsNear(3.1e-8_(1) * Metre)),
                AbsoluteErrorFrom(0.2 * Metre, IsNear(2.3e-8_(1) * Metre))),
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

    EXPECT_EQ(1372, function_invocations);
    EXPECT_EQ(470, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                AbsoluteErrorFrom(0.114589 * Metre, IsNear(3.4e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.555649 * Metre, IsNear(2.7e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.852547 * Metre,
                                  IsNear(3.9e-7_(1) * Metre))),
            Componentwise(
                AbsoluteErrorFrom(0.109337 * Metre, IsNear(7.0e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.860556 * Metre, IsNear(4.7e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.564135 * Metre,
                                  IsNear(3.4e-7_(1) * Metre))),
            Componentwise(
                AbsoluteErrorFrom(0.688823 * Metre, IsNear(2.4e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.117274 * Metre, IsNear(5.8e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.267465 * Metre,
                                  IsNear(1.3e-6_(1) * Metre)))));
  }
  function_invocations = 0;
  gradient_invocations = 0;
  {
    auto const minima =
        optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                   /*number_of_rounds=*/std::nullopt,
                                   tolerance);

    EXPECT_EQ(199, function_invocations);
    EXPECT_EQ(124, gradient_invocations);
    EXPECT_THAT(
        minima,
        ElementsAre(
            Componentwise(
                AbsoluteErrorFrom(0.688823 * Metre, IsNear(2.6e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.117274 * Metre, IsNear(4.1e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.267465 * Metre,
                                  IsNear(1.2e-6_(1) * Metre))),
            Componentwise(
                AbsoluteErrorFrom(0.114589 * Metre, IsNear(3.9e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.555649 * Metre, IsNear(6.3e-7_(1) * Metre)),
                AbsoluteErrorFrom(0.852547 * Metre,
                                  IsNear(4.8e-7_(1) * Metre)))));
  }
}

// A function that looks like the opposite of the gravitational potential.
TEST_F(GlobalOptimizationTest, Potential) {
  mathematica::Logger logger(TEMP_DIR / "potential.wl");
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
    logger.Append("minima", minima, mathematica::ExpressIn(Metre));
    logger.Flush();

    EXPECT_EQ(29638, function_invocations);
    EXPECT_EQ(28689, gradient_invocations);
    EXPECT_THAT(minima, IsEmpty());
  }
  function_invocations = 0;
  gradient_invocations = 0;
  {
    auto const minima =
        optimizer.FindGlobalMinima(/*points_per_round=*/10,
                                   /*number_of_rounds=*/std::nullopt,
                                   tolerance);

    EXPECT_EQ(5101, function_invocations);
    EXPECT_EQ(5065, gradient_invocations);
    EXPECT_THAT(minima, IsEmpty());
  }
}

}  // namespace numerics
}  // namespace principia
