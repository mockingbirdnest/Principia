#include "numerics/gradient_descent.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {

using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::Vector;
using quantities::Exponentiation;
using quantities::Inverse;
using quantities::Length;
using quantities::Pow;
using quantities::si::Micro;
using quantities::si::Metre;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::IsNear;
using testing_utilities::operator""_;

class GradientDescentTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
};

TEST_F(GradientDescentTest, Quadratic) {
  auto field = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return Pow<2>(coordinates.x - 1 * Metre) +
           Pow<2>(coordinates.y - 2 * Metre) +
           Pow<2>(coordinates.z + 3 * Metre);
  };
  auto gradient = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return Vector<Length, World>(
        {2 * (coordinates.x - 1 * Metre),
         2 * (coordinates.y - 2 * Metre),
         2 * (coordinates.z + 3 * Metre)});
  };

  Position<World> const expected_minimum =
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre});
  auto const actual_minimum =
      BroydenFletcherGoldfarbShanno<Exponentiation<Length, 2>, Position<World>>(
          /*start_argument=*/World::origin,
          field,
          gradient,
          /*tolerance=*/1 * Micro(Metre));
  EXPECT_THAT(actual_minimum, AlmostEquals(expected_minimum, 2));
}

TEST_F(GradientDescentTest, Quartic) {
  auto field = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return Pow<4>(coordinates.x - 1 * Metre) +
           Pow<4>(coordinates.y - 2 * Metre) +
           Pow<4>(coordinates.z + 3 * Metre);
  };
  auto gradient = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return Vector<Exponentiation<Length, 3>, World>(
        {4 * Pow<3>(coordinates.x - 1 * Metre),
         4 * Pow<3>(coordinates.y - 2 * Metre),
         4 * Pow<3>(coordinates.z + 3 * Metre)});
  };

  Position<World> const expected_minimum =
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre});
  auto const actual_minimum =
      BroydenFletcherGoldfarbShanno<Exponentiation<Length, 4>, Position<World>>(
          /*start_argument=*/World::origin,
          field,
          gradient,
          /*tolerance=*/1 * Micro(Metre));
  EXPECT_THAT(
      actual_minimum,
      AbsoluteErrorFrom(expected_minimum, IsNear(437_(1) * Micro(Metre))));
}

TEST_F(GradientDescentTest, Gaussian) {
  static constexpr auto ω = 1 / Metre;
  auto field = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return -std::exp(-Pow<2>(ω * (coordinates.x - 1 * Metre))) -
           std::exp(-Pow<2>(ω * (coordinates.y - 2 * Metre))) -
           std::exp(-Pow<2>(ω * (coordinates.z + 3 * Metre)));
  };
  auto gradient = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return Vector<Inverse<Length>, World>(
        {2 * Pow<2>(ω) * (coordinates.x - 1 * Metre) *
             std::exp(-Pow<2>(ω * (coordinates.x - 1 * Metre))),
         2 * Pow<2>(ω) * (coordinates.y - 2 * Metre) *
             std::exp(-Pow<2>(ω * (coordinates.y - 2 * Metre))),
         2 * Pow<2>(ω) * (coordinates.z + 3 * Metre) *
             std::exp(-Pow<2>(ω * (coordinates.z + 3 * Metre)))});
  };

  Position<World> const expected_minimum =
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre});
  auto const actual_minimum =
      BroydenFletcherGoldfarbShanno<double, Position<World>>(
          /*start_argument=*/World::origin,
          field,
          gradient,
          /*tolerance=*/1 * Micro(Metre));
  EXPECT_THAT(
      actual_minimum,
      AbsoluteErrorFrom(expected_minimum, IsNear(0.82_(1) * Micro(Metre))));
}

TEST_F(GradientDescentTest, Rosenbrock) {
  // See [NW06] exercise 3.1 for the definition of the function.
  static constexpr auto scale = 1 * Metre;
  auto field = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates() / scale;
    return 100 * Pow<2>(coordinates.y - Pow<2>(coordinates.x)) +
           Pow<2>(1 - coordinates.x);
  };
  auto gradient = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates() / scale;
    return Vector<Inverse<Length>, World>(
        {(-400 * coordinates.x * (coordinates.y - Pow<2>(coordinates.x)) -
          2 * (1 - coordinates.x)) /
             scale,
         (200 * (coordinates.y - Pow<2>(coordinates.x))) / scale,
         Inverse<Length>{}});
  };

  // See [NW06] exercise 3.1 for the starting positions.
  Position<World> const expected_minimum =
      World::origin + Displacement<World>({scale, scale, 0 * Metre});
  {
    auto const actual_minimum =
        BroydenFletcherGoldfarbShanno<double, Position<World>>(
            /*start_argument=*/World::origin +
                Displacement<World>({1.2 * scale, 1.2 * scale, 0 * Metre}),
            field,
            gradient,
            /*tolerance=*/1 * Micro(Metre));
    EXPECT_THAT(
        actual_minimum,
        AbsoluteErrorFrom(expected_minimum, IsNear(0.96_(1) * Micro(Metre))));
  }
  {
    auto const actual_minimum =
        BroydenFletcherGoldfarbShanno<double, Position<World>>(
            /*start_argument=*/World::origin +
                Displacement<World>({-1.2 * scale, scale, 0 * Metre}),
            field,
            gradient,
            /*tolerance=*/1 * Micro(Metre));
    EXPECT_THAT(
        actual_minimum,
        AbsoluteErrorFrom(expected_minimum, IsNear(0.047_(1) * Micro(Metre))));
  }
}

}  // namespace numerics
}  // namespace principia
