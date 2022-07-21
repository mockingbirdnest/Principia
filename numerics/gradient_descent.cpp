#include "numerics/gradient_descent.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Position;
using quantities::Exponentiation;
using quantities::Pow;
using quantities::si::Metre;

class GradientDescentTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;
};

TEST_F(GradientDescentTest, Basic) {
  auto field = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return Pow<2>(coordinates.x - 1 * Metre) *
           Pow<4>(coordinates.y - 2 * Metre) *
           Pow<6>(coordinates.z + 3 * Metre);
  };
  auto gradient = [](Position<World> const& position) {
    auto const coordinates = (position - World::origin).coordinates();
    return 2 * (coordinates.x - 1 * Metre) *
           Pow<4>(coordinates.y - 2 * Metre) *
           Pow<6>(coordinates.z + 3 * Metre) +
           Pow<2>(coordinates.x - 1 * Metre) *
           4 * (coordinates.y - 2 * Metre) *
           Pow<6>(coordinates.z + 3 * Metre) +
           Pow<2>(coordinates.x - 1 * Metre) *
           Pow<4>(coordinates.y - 2 * Metre) *
           6 * (coordinates.z + 3 * Metre);
  };

  auto const minimum = GradientDescent<Exponentiation<Length, 12>>(
      /*start_position=*/World::origin, field, gradient, []() { return true; });
}

}  // namespace numerics
}  // namespace principia
