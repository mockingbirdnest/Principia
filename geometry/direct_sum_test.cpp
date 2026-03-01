#include "geometry/direct_sum.hpp"

#include <chrono>
#include <tuple>

#include "base/algebra.hpp"
#include "base/tags.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/point.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/algebra.hpp"

namespace principia {
namespace geometry {

using namespace principia::base::_algebra;
using namespace principia::base::_tags;
using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_point;
using namespace principia::geometry::_space;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_algebra;

using ℝ² = DirectSum<double, double>;
using World = Frame<serialization::Frame::TestTag,
                    Inertial,
                    Handedness::Right,
                    serialization::Frame::TEST>;

// A helper type for counting constructions.  It must be an affine type.
struct ConstructionCounter {
  static std::int64_t initialized_count;
  static std::int64_t uninitialized_count;

  constexpr explicit ConstructionCounter(uninitialized_t) {
    ++uninitialized_count;
  }

  constexpr ConstructionCounter() {
    ++initialized_count;
  }

  ConstructionCounter& operator+=(double) {
    return *this;
  };
  ConstructionCounter& operator-=(double) {
    return *this;
  };
};

std::int64_t ConstructionCounter::initialized_count = 0;
std::int64_t ConstructionCounter::uninitialized_count = 0;

double operator-(ConstructionCounter const& left,
                 ConstructionCounter const& right) {
  return 0.0;
}

ConstructionCounter operator+(ConstructionCounter const& left, double right) {
  return left;
}

ConstructionCounter operator-(ConstructionCounter const& left, double right) {
  return left;
}

static_assert(affine<ConstructionCounter>);


TEST(DirectSumTest, AlgebraConcepts) {
  static_assert(affine<DirectSum<std::byte*, double>>);

  static_assert(additive_group<DirectSum<std::chrono::seconds, double>>);
  static_assert(
      homogeneous_module<DirectSum<std::chrono::seconds, double>, int>);

  static_assert(real_affine_space<DirectSum<Point<Length>>>);
  static_assert(real_affine_space<DirectSum<Point<Length>, Speed>>);

  static_assert(real_vector_space<ℝ²>);
  static_assert(real_vector_space<DirectSum<Length, Mass, Time>>);

  static_assert(homogeneous_vector_space<DirectSum<Length>, Time>);

  static_assert(vector_space<DirectSum<IntegerModulo<2>, IntegerModulo<2>>,
                             IntegerModulo<2>>);

  static_assert(hilbert<ℝ²>);
  static_assert(!hilbert<DirectSum<IntegerModulo<2>, IntegerModulo<2>>>);
}

// Helper concepts for Unordered test.
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

TEST(DirectSumTest, Unordered) {
  static_assert(!less_than<ℝ²>);
  static_assert(!greater_than<ℝ²>);
  static_assert(!less_than_or_equal<ℝ²>);
  static_assert(!greater_than_or_equal<ℝ²>);
}

// Helper concepts for DegreesOfFreedomIsMerelyAffine test.
template<typename T>
concept has_unary_plus = requires(T a) { +a; };
template<typename T>
concept has_unary_minus = requires(T a) { -a; };
template<typename T, typename U>
concept has_plus = requires(T a, U b) { a + b; };
template<typename T, typename U>
concept has_minus = requires(T a, U b) { a - b; };
template<typename T, typename U>
concept has_times = requires(T a, U b) { a * b; };
template<typename T, typename U>
concept has_divided_by = requires(T a, U b) { a / b; };

TEST(DirectSumTest, DegreesOfFreedomIsMerelyAffine) {
  using DegreesOfFreedom = DirectSum<Position<World>, Velocity<World>>;
  static_assert(real_affine_space<DegreesOfFreedom>);

  static_assert(!has_unary_plus<DegreesOfFreedom>);
  static_assert(!has_unary_minus<DegreesOfFreedom>);
  static_assert(!has_plus<DegreesOfFreedom, DegreesOfFreedom>);
  static_assert(!has_minus<Difference<DegreesOfFreedom>, DegreesOfFreedom>);
  static_assert(!has_times<DegreesOfFreedom, double>);
  static_assert(!has_times<double, DegreesOfFreedom>);
  static_assert(!has_divided_by<DegreesOfFreedom, double>);
  static_assert(!has_divided_by<double, DegreesOfFreedom>);
}

TEST(DirectSumTest, Constructors) {
  EXPECT_EQ(DirectSum<double>(), DirectSum{0.0});
  EXPECT_EQ(DirectSum<Length>(), DirectSum{0 * Metre});
}

TEST(DirectSumTest, UninitializedConstruction) {
  DirectSum<double, Length, ConstructionCounter> d1(uninitialized);
  EXPECT_EQ(get<1>(d1), 0 * Metre);
  EXPECT_EQ(0, ConstructionCounter::initialized_count);
  EXPECT_EQ(1, ConstructionCounter::uninitialized_count);
  DirectSum<double, Length, ConstructionCounter> d2;
  EXPECT_EQ(get<1>(d2), 0 * Metre);
  EXPECT_EQ(1, ConstructionCounter::initialized_count);
  EXPECT_EQ(1, ConstructionCounter::uninitialized_count);
}

TEST(DirectSumTest, UnaryPlus) {
  EXPECT_EQ(+ℝ²(1, 2), ℝ²(1, 2));
}

TEST(DirectSumTest, Negation) {
  EXPECT_EQ(-ℝ²(1, 2), ℝ²(-1, -2));
}

TEST(DirectSumTest, Addition) {
  EXPECT_EQ(ℝ²(1, 2) + ℝ²(3, 4), ℝ²(4, 6));

  // Affine addition.
  EXPECT_EQ(DirectSum<Length>(1 * Metre) + DirectSum<Point<Length>>(),
            DirectSum<Point<Length>>(Point<Length>() + 1 * Metre));
  EXPECT_EQ(DirectSum<Point<Length>>() + DirectSum<Length>(1 * Metre),
            DirectSum<Point<Length>>(Point<Length>() + 1 * Metre));
}

TEST(DirectSumTest, Subtraction) {
  EXPECT_EQ(ℝ²(1, 2) - ℝ²(3, 4), ℝ²(-2, -2));
}

TEST(DirectSumTest, Multiplication) {
  EXPECT_EQ(2 * ℝ²(3, 4), ℝ²(6, 8));
  EXPECT_EQ(ℝ²(3, 4) * 2, ℝ²(6, 8));
}

TEST(DirectSumTest, Division) {
  // Note: ℝ²(3, 4) / 2 doesn't work because integers are not a field.
  EXPECT_EQ(ℝ²(3, 4) / 2.0, ℝ²(1.5, 2));
}

TEST(DirectSumTest, InnerProduct) {
  ℝ² one_two = {1, 2};
  ℝ² three_four = {3, 4};

  EXPECT_EQ(InnerProduct(one_two, three_four), 11);
  EXPECT_EQ(three_four.Norm²(), 25);
  EXPECT_EQ(three_four.Norm(), 5);
}

TEST(DirectSumTest, InnerProductOfDirectSumOfVector) {
  DirectSum<Vector<double, World>> one_two_three = {
      Vector<double, World>({1, 2, 3})};
  DirectSum<Vector<double, World>> four_five_six = {
      Vector<double, World>({4, 5, 6})};

  EXPECT_EQ(InnerProduct(one_two_three, four_five_six), 32);
}

TEST(DirectSumTest, Get) {
  ℝ² one_two = {1, 2};
  EXPECT_EQ(get<0>(one_two), 1);
}

TEST(DirectSumTest, StructuredBindingsConst) {
  DirectSum<Length, Time> one_two = {1 * Metre, 2 * Second};

  auto const& [length, time] = one_two;
  EXPECT_EQ(length, 1 * Metre);
  EXPECT_EQ(time, 2 * Second);
}

TEST(DirectSumTest, StructuredBindingsNonConst) {
  DirectSum<Length, Time> one_two = {1 * Metre, 2 * Second};

  auto& [length, time] = one_two;
  length += 1 * Metre;
  time += 4 * Second;

  EXPECT_EQ(one_two, (DirectSum<Length, Time>{2 * Metre, 6 * Second}));
}

TEST(DirectSumTest, FixedVector) {
  FixedVector<ℝ², 1>(std::array<ℝ², 1>{ℝ²{1, 2}});
}

}  // namespace geometry
}  // namespace principia
