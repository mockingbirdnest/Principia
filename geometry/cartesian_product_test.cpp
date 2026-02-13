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

TEST(CartesianProductTest, AlgebraConcepts) {
  static_assert(affine<DirectSum<std::byte*, double>>);

  static_assert(additive_group<DirectSum<std::chrono::seconds, double>>);
  static_assert(
      homogeneous_module<DirectSum<std::chrono::seconds, double>, int>);

  static_assert(real_affine_space<DirectSum<Point<Length>>>);
  static_assert(real_affine_space<DirectSum<Point<Length>, Speed>>);

  static_assert(real_vector_space<DirectSum<double, double>>);
  static_assert(real_vector_space<DirectSum<Length, Mass, Time>>);

  using ℝ² = DirectSum<double, double>;
  static_assert(hilbert<ℝ², ℝ²>);
}

TEST(CartesianProductTest, Constructors) {
  EXPECT_EQ(DirectSum<double>(), DirectSum{0.0});
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

TEST(CartesianProductTest, InnerProduct) {
  DirectSum<double, double> one_two = {1, 2};
  DirectSum<double, double> three_four = {3, 4};

  EXPECT_EQ(InnerProduct(one_two, three_four), 11);
  EXPECT_EQ(three_four.Norm²(), 25);
  EXPECT_EQ(three_four.Norm(), 5);
}

TEST(CartesianProductTest, FixedVector) {
  FixedVector<DirectSum<double, double>, 1>(
      std::array<DirectSum<double, double>, 1>{
          DirectSum<double, double>{1, 2}});
}

}  // namespace geometry
}  // namespace principia
