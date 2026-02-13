#include "geometry/cartesian_product.hpp"

#include <tuple>

#include "base/algebra.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace geometry {

using namespace principia::base::_algebra;
using namespace principia::geometry::_cartesian_product;
using namespace principia::geometry::_point;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

TEST(CartesianProductTest, Concepts) {
  static_assert(real_vector_space<DirectSum<double>>);
  static_assert(real_vector_space<DirectSum<double, double>>);
  static_assert(real_vector_space<DirectSum<Length>>);
  static_assert(real_vector_space<DirectSum<Length, Mass, Time>>);

  static_assert(real_affine_space<DirectSum<Point<Length>>>);
  static_assert(real_affine_space<DirectSum<Point<Length>, Speed>>);
}

TEST(CartesianProductTest, Constructors) {
  EXPECT_EQ(DirectSum<Length>(), DirectSum{0 * Metre});
  EXPECT_EQ(DirectSum<Length>(std::tuple<Length>(4 * Metre)),
            DirectSum{4 * Metre});
}

TEST(CartesianProductTest, StructuredBindingsConst) {
  DirectSum<Length, Time> one_two = {1 * Metre, 2 * Second};

  auto const& [length, time] = one_two;
  EXPECT_EQ(length, 1 * Metre);
  EXPECT_EQ(time, 2 * Second);
}

TEST(CartesianProductTest, StructuredBindingsNonConst) {
  DirectSum<Length, Time> one_two = {1 * Metre, 2 * Second};

  auto& [length, time] = one_two;
  length += 1 * Metre;
  time += 4 * Second;

  EXPECT_EQ(one_two, (DirectSum<Length, Time>{2 * Metre, 6 * Second}));
}

TEST(CartesianProductTest, FixedVector) {
  FixedVector<DirectSum<double, double>, 1>(
      std::array<DirectSum<double, double>, 1>{
          DirectSum<double, double>{1, 2}});
}

}  // namespace geometry
}  // namespace principia
