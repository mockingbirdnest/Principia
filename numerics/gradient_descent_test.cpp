#include "numerics/gradient_descent.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For IsOkAndHolds.
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {

using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_space;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_gradient_descent;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

class GradientDescentTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag>;
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
  auto directional_gradient = [](Position<World> const& position,
                                 Displacement<World> const& displacement) {
    auto const position_coordinates = (position - World::origin).coordinates();
    auto const displacement_coordinates = displacement.coordinates();
    return 2 *
           ((position_coordinates.x - 1 * Metre) * displacement_coordinates.x +
            (position_coordinates.y - 2 * Metre) * displacement_coordinates.y +
            (position_coordinates.z + 3 * Metre) * displacement_coordinates.z);
  };

  Position<World> const expected_minimum =
      World::origin + Displacement<World>({1 * Metre, 2 * Metre, -3 * Metre});
  {
    auto const actual_minimum =
        BroydenFletcherGoldfarbShanno<Exponentiation<Length, 2>,
                                      Position<World>>(
            /*start_argument=*/World::origin,
            field,
            gradient,
            /*tolerance=*/1 * Micro(Metre));
    EXPECT_THAT(actual_minimum,
                IsOkAndHolds(AlmostEquals(expected_minimum, 2)));
  }

  {
    auto const actual_minimum =
        BroydenFletcherGoldfarbShanno<Exponentiation<Length, 2>,
                                      Position<World>>(
            /*start_argument=*/World::origin,
            field,
            gradient,
            directional_gradient,
            /*tolerance=*/1 * Micro(Metre));
    EXPECT_THAT(actual_minimum,
                IsOkAndHolds(AlmostEquals(expected_minimum, 2)));
  }
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
  EXPECT_THAT(actual_minimum,
              IsOkAndHolds(AbsoluteErrorFrom(expected_minimum,
                                             IsNear(437_(1) * Micro(Metre)))));
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
  EXPECT_THAT(actual_minimum,
              IsOkAndHolds(AbsoluteErrorFrom(expected_minimum,
                                             IsNear(0.78_(1) * Micro(Metre)))));
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
    EXPECT_THAT(actual_minimum,
                IsOkAndHolds(AbsoluteErrorFrom(
                    expected_minimum, IsNear(0.69_(1) * Micro(Metre)))));
  }
  {
    auto const actual_minimum =
        BroydenFletcherGoldfarbShanno<double, Position<World>>(
            /*start_argument=*/World::origin +
                Displacement<World>({-1.2 * scale, scale, 0 * Metre}),
            field,
            gradient,
            /*tolerance=*/1 * Micro(Metre));
    EXPECT_THAT(actual_minimum,
                IsOkAndHolds(AbsoluteErrorFrom(
                    expected_minimum, IsNear(0.86_(1) * Micro(Metre)))));
  }
}

TEST_F(GradientDescentTest, Fixed) {
  using Scalar = Voltage;
  using Argument = FixedVector<Length, 2>;

  auto field = [](Argument const& argument) -> Voltage {
    return 3 * (Volt / Pow<2>(Metre)) * (Pow<2>(argument[0] - 1 * Metre) +
                                         Pow<2>(argument[1] + 2 * Metre));
  };
  auto gradient = [](Argument const& argument) {
    return FixedVector<Quotient<Voltage, Length>, 2>(
        {6 * (Volt / Pow<2>(Metre)) * (argument[0] - 1 * Metre),
         6 * (Volt / Pow<2>(Metre)) * (argument[1] + 2 * Metre)});
  };

  Argument const expected_minimum({1 * Metre, -2 * Metre});
  auto const actual_minimum =
      BroydenFletcherGoldfarbShanno<Scalar, Argument>(
          /*start_argument=*/Argument({0 * Metre, 0 * Metre}),
          field,
          gradient,
          /*tolerance=*/1 * Micro(Metre));
  EXPECT_THAT(actual_minimum, IsOkAndHolds(AlmostEquals(expected_minimum, 1)));
}

}  // namespace numerics
}  // namespace principia
