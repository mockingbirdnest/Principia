#include "geometry/cartesian_product.hpp"

#include <chrono>
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

using ℝ² = DirectSum<double, double>;

TEST(CartesianProductTest, AlgebraConcepts) {
  static_assert(affine<DirectSum<std::byte*, double>>);

  static_assert(additive_group<DirectSum<std::chrono::seconds, double>>);
  static_assert(
      homogeneous_module<DirectSum<std::chrono::seconds, double>, int>);

  static_assert(real_affine_space<DirectSum<Point<Length>>>);
  static_assert(real_affine_space<DirectSum<Point<Length>, Speed>>);

  static_assert(real_vector_space<DirectSum<double, double>>);
  static_assert(real_vector_space<DirectSum<Length, Mass, Time>>);

  static_assert(homogeneous_vector_space<DirectSum<Length>, Time>);

  static_assert(hilbert<ℝ², ℝ²>);
}

template<typename T>
concept less_than = requires(T a, T b) {
  { a < b };
};
template<typename T>
concept greater_than = requires(T a, T b) {
  { a > b };
};
template<typename T>
concept less_than_or_equal = requires(T a, T b) {
  { a <= b };
};
template<typename T>
concept greater_than_or_equal = requires(T a, T b) {
  { a >= b };
};

TEST(CartesianProductTest, Unordered) {
  static_assert(!less_than<ℝ²>);
}

TEST(CartesianProductTest, Constructors) {
  EXPECT_EQ(DirectSum<double>(), DirectSum{0.0});
  EXPECT_EQ(DirectSum<Length>(), DirectSum{0 * Metre});
  EXPECT_EQ(DirectSum<Length>(std::tuple<Length>(4 * Metre)),
            DirectSum{4 * Metre});
}

TEST(CartesianProductTest, UnaryPlus) {
  EXPECT_EQ(+ℝ²(1, 2), ℝ²(1, 2));
}

TEST(CartesianProductTest, Negation) {
  EXPECT_EQ(-ℝ²(1, 2), ℝ²(-1, -2));
}

TEST(CartesianProductTest, Addition) {
  EXPECT_EQ(ℝ²(1, 2) + ℝ²(3, 4), ℝ²(4, 6));
}

TEST(CartesianProductTest, Subtraction) {
  EXPECT_EQ(ℝ²(1, 2) - ℝ²(3, 4), ℝ²(-2, -2));
}

TEST(CartesianProductTest, Multiplication) {
  EXPECT_EQ(2 * ℝ²(3, 4), ℝ²(6, 8));
  EXPECT_EQ(ℝ²(3, 4) * 2, ℝ²(6, 8));
}

TEST(CartesianProductTest, Division) {
  // Note: ℝ²(3, 4) / 2 doesn't work because integers are not a field.
  EXPECT_EQ(ℝ²(3, 4) / 2.0, ℝ²(1.5, 2));
}

TEST(CartesianProductTest, InnerProduct) {
  DirectSum<double, double> one_two = {1, 2};
  DirectSum<double, double> three_four = {3, 4};

  EXPECT_EQ(InnerProduct(one_two, three_four), 11);
  EXPECT_EQ(three_four.Norm²(), 25);
  EXPECT_EQ(three_four.Norm(), 5);
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
