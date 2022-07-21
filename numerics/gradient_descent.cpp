#include "numerics/gradient_descent.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace numerics {

using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::Vector;
using quantities::Length;
using quantities::Pow;
using quantities::Exponentiation;
using quantities::si::Micro;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;

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
  auto const actual_minimum = GradientDescent<Exponentiation<Length, 2>, World>(
      /*start_position=*/World::origin,
      field,
      gradient,
      /*tolerance=*/1 * Micro(Metre));
  EXPECT_THAT(actual_minimum, AlmostEquals(expected_minimum, 0));
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
  auto const actual_minimum = GradientDescent<Exponentiation<Length, 4>, World>(
      /*start_position=*/World::origin,
      field,
      gradient,
      /*tolerance=*/1 * Micro(Metre));
  EXPECT_THAT(actual_minimum, AlmostEquals(expected_minimum, 2));
}

}  // namespace numerics
}  // namespace principia
