#pragma once

#include "gtest/gtest.h"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
void TestEquality(T const& low, T const& high) {
  EXPECT_TRUE(low == low) << "low = " << low;
  EXPECT_TRUE(high == high) << "high = " << high;
  EXPECT_TRUE(high != low) << "high = " << high << ", low = " << low;
  EXPECT_TRUE(low != high) << "low = " << low << ", high = " << high;

  EXPECT_FALSE(high == low) << "high = " << high << ", low = " << low;
  EXPECT_FALSE(low == high) << "low = " << low << ", high = " << high;
  EXPECT_FALSE(low != low) << "low = " << low;
  EXPECT_FALSE(high != high) << "high = " << high;
}

template<typename T>
void TestOrder(T const& low, T const& high) {
  TestEquality(low, high);
  EXPECT_TRUE(high > low) << "high = " << high << ", low = " << low;
  EXPECT_TRUE(low < high) << "low = " << low << ", high = " << high;
  EXPECT_TRUE(low >= low) << "low = " << low;
  EXPECT_TRUE(low <= low) << "low = " << low;
  EXPECT_TRUE(high >= high) << "high = " << high;
  EXPECT_TRUE(high <= high) << "high = " << high;
  EXPECT_TRUE(high >= low) << "high = " << high << ", low = " << low;
  EXPECT_TRUE(low <= high) << "low = " << low << ", high = " << high;

  EXPECT_FALSE(low > low) << "low = " << low;
  EXPECT_FALSE(low < low) << "low = " << low;
  EXPECT_FALSE(high > high) << "high = " << high;
  EXPECT_FALSE(high < high) << "high = " << high;
  EXPECT_FALSE(low > high) << "low = " << low << ", high = " << high;
  EXPECT_FALSE(high < low) << "high = " << high << ", low = " << low;
  EXPECT_FALSE(low >= high) << "low = " << low << ", high = " << high;
  EXPECT_FALSE(high <= low) << "high = " << high << ", low = " << low;
}

template<typename T>
void TestAdditiveGroup(T const& zero, T const& a, T const& b, T const& c,
                       std::int64_t const max_ulps) {
  EXPECT_EQ(a, +a);
  EXPECT_EQ(a + zero, a);
  EXPECT_EQ(zero + b, b);
  EXPECT_EQ(a - a, zero);
  EXPECT_EQ(-a - b, -(a + b));
  EXPECT_THAT((a + b) + c, AlmostEquals(a + (b + c), max_ulps));
  EXPECT_THAT(a - b - c, AlmostEquals(a - (b + c), max_ulps));
  EXPECT_EQ(a + b, b + a);
  T accumulator = zero;
  accumulator += a;
  accumulator += b;
  accumulator -= c;
  EXPECT_THAT(accumulator, AlmostEquals(a + b - c, max_ulps));
}


template<typename T>
void TestGroup(T const& identity, T const& a, T const& b, T const& c,
               T (*operation)(T const&, T const&), T (*inverse)(T const&),
               std::int64_t const max_ulps) {
  EXPECT_EQ(operation(a, identity), a);
  EXPECT_EQ(operation(b, identity), b);
  EXPECT_THAT(operation(a, inverse(a)), AlmostEquals(identity, max_ulps));
  EXPECT_THAT(operation(inverse(b), inverse(a)),
              AlmostEquals(inverse(operation(a, b)), max_ulps));
  EXPECT_THAT(operation(operation(a, b), c),
              AlmostEquals(operation(a, operation(b, c)), max_ulps));
  EXPECT_THAT(operation(operation(a, inverse(b)), inverse(c)),
              AlmostEquals(operation(a, inverse(operation(c, b))), max_ulps));
}

template<typename T>
void TestAbelianMultiplicativeGroup(T const& one, T const& a, T const& b, T const& c,
                             std::int64_t const max_ulps) {
  EXPECT_EQ(a * one, a);
  EXPECT_EQ(one * b, b);
  EXPECT_EQ(a / a, one);
  EXPECT_THAT((one / a) / b, AlmostEquals(one / (a * b), max_ulps));
  EXPECT_THAT((a * b) * c, AlmostEquals(a * (b * c), max_ulps));
  EXPECT_THAT(a / b / c, AlmostEquals(a / (b * c), max_ulps));
  EXPECT_EQ(a * b, b * a);
  T accumulator = one;
  accumulator *= a;
  accumulator *= b;
  accumulator /= c;
  EXPECT_THAT(accumulator, AlmostEquals(a * b / c, max_ulps));
}

template<typename T>
void TestNonAbelianMultiplicativeGroup(T const& one, T const& a, T const& b, 
                                       T const& c, std::int64_t const max_ulps) {
  EXPECT_EQ(a * one, a);
  EXPECT_EQ(one * b, b);
  EXPECT_EQ(a / a, one);
  EXPECT_THAT((one / a) / b, AlmostEquals(one / (a * b), max_ulps));
  EXPECT_THAT((a * b) * c, AlmostEquals(a * (b * c), max_ulps));
  EXPECT_THAT(a / b / c, AlmostEquals(a / (c * b), max_ulps));
  T accumulator = one;
  accumulator *= a;
  accumulator *= b;
  accumulator /= c;
  EXPECT_THAT(accumulator, AlmostEquals(a * b / c, max_ulps));
}

// The Greek letters cause a warning when stringified by the macros, because
// apparently Visual Studio doesn't encode strings in UTF-8 by default.
#pragma warning(disable: 4566)

template<typename Map, typename Scalar, typename U, typename V>
void TestBilinearMap(Map const& map, U const& u1, U const& u2, V const& v1,
                     V const& v2, Scalar const& λ,
                     std::int64_t const max_ulps) {
  EXPECT_THAT(map(u1 + u2, v1),
              AlmostEquals(map(u1, v1) + map(u2, v1), max_ulps));
  EXPECT_THAT(map(u1, v1 + v2),
              AlmostEquals(map(u1, v1) + map(u1, v2), max_ulps));
  EXPECT_THAT(map(λ * u1, v1), AlmostEquals(map(u1, λ * v1), max_ulps));
  EXPECT_THAT(λ * map(u1, v1), AlmostEquals(map(u1, λ * v1), max_ulps));
  EXPECT_THAT(map(u2 * λ, v2), AlmostEquals(map(u2, v2 * λ), max_ulps));
  EXPECT_THAT(map(u2, v2) * λ, AlmostEquals(map(u2 * λ, v2), max_ulps));
}

template<typename Map, typename Scalar, typename U>
void TestSymmetricBilinearMap(Map const& map, U const& u1, U const& u2,
                               U const& v1, U const& v2, Scalar const& λ,
                               std::int64_t const max_ulps) {
  TestBilinearMap(map, u1, u2, v1, v2, λ, max_ulps);
  EXPECT_THAT(map(u1, v1), AlmostEquals(map(v1, u1), max_ulps));
  EXPECT_THAT(map(u2, v2), AlmostEquals(map(v2, u2), max_ulps));
}

template<typename Map, typename Scalar, typename U>
void TestSymmetricPositiveDefiniteBilinearMap(
    Map const& map,
    U const& u1,
    U const& u2, U const& v1,
    U const& v2, Scalar const& λ,
    std::int64_t const max_ulps) {
  TestSymmetricBilinearMap(map, u1, u2, v1, v2, λ, max_ulps);
  auto zero = map(u1, u1) - map(u1, u1);
  EXPECT_TRUE(map(u1, u1) > zero) << map(u1, u1);
  EXPECT_TRUE(map(u2, u2) > zero) << map(u2, u2);
  EXPECT_TRUE(map(v1, v1) > zero) << map(v1, v1);
  EXPECT_TRUE(map(v2, v2) > zero) << map(v2, v2);
}

template<typename Map, typename Scalar, typename U>
void TestAlternatingBilinearMap(Map const& map, U const& u1, U const& u2,
                                U const& v1, U const& v2, Scalar const& λ,
                                std::int64_t const max_ulps) {
  TestBilinearMap(map, u1, u2, v1, v2, λ, max_ulps);
  auto zero = map(u1, u1) - map(u1, u1);
  EXPECT_THAT(map(u1, u1), AlmostEquals(zero, max_ulps));
  EXPECT_THAT(map(u2, u2), AlmostEquals(zero, max_ulps));
  EXPECT_THAT(map(v1, v2), AlmostEquals(-map(v2, v1), max_ulps));
}

template<typename Map, typename Scalar, typename U>
void TestLieBracket(Map const& map, U const& u1, U const& u2, U const& v1,
                    U const& v2, Scalar const& λ,
                    std::int64_t const max_ulps) {
  TestAlternatingBilinearMap(map, u1, u2, v1, v2, λ, max_ulps);
  EXPECT_THAT(map(u1, map(u2, v1)) + map(u2, map(v1, u1)),
              AlmostEquals(-map(v1, map(u1, u2)), max_ulps));
}

template<typename Vector, typename Scalar>
void TestVectorSpace(Vector const& nullVector, Vector const& u, Vector const& v,
                     Vector const& w, Scalar const& zero, Scalar const& unit,
                     Scalar const& α, Scalar const& β,
                     std::int64_t const max_ulps) {
  TestAdditiveGroup(nullVector, u, v, w, max_ulps);
  EXPECT_THAT((α * β) * v, AlmostEquals(α * (β * v), max_ulps));
  EXPECT_EQ(unit * w, w);
  EXPECT_EQ(u / unit, u);
  EXPECT_EQ(zero * u, nullVector);
  EXPECT_THAT(β * (u + v), AlmostEquals(β * u + v * β, max_ulps));
  EXPECT_THAT((w + v) / α, AlmostEquals(w / α + v / α, max_ulps));
  EXPECT_THAT((α + β) * w, AlmostEquals(α * w + β * w, max_ulps));
  EXPECT_THAT(v * (α + β), AlmostEquals(α * v + β * v, max_ulps));
  Vector vector = u;
  vector *= α;
  EXPECT_EQ(α * u, vector);
  vector /= α;
  EXPECT_THAT(u, AlmostEquals(vector, max_ulps));
  vector *= zero;
  EXPECT_EQ(vector, nullVector);
}

template<typename Vector, typename Scalar, typename Map>
void TestInnerProductSpace(Map const& map, Vector const& nullVector,
                           Vector const& u, Vector const& v, Vector const& w,
                           Vector const& a, Scalar const& zero,
                           Scalar const& unit, Scalar const& α, Scalar const& β,
                           std::int64_t const max_ulps) {
  TestVectorSpace(nullVector, u, v, w, zero, unit, α, β, max_ulps);
  TestSymmetricPositiveDefiniteBilinearMap(map, u, v, w, a, α, max_ulps);
}

#pragma warning(default: 4566)

template<typename T>
void TestField(T const& zero, T const& one, T const& a, T const& b, T const& c,
               T const& x, T const& y, std::int64_t const max_ulps) {
  TestAdditiveGroup(zero, a, b, c, max_ulps);
  TestAbelianMultiplicativeGroup(one, c, x, y, max_ulps);
  TestVectorSpace(zero, a, b, c, zero, one, x, y, max_ulps);
}

template<typename T>
void TestSkewField(T const& zero, T const& one, T const& a, T const& b, T const& c,
                   T const& x, T const& y, std::int64_t const max_ulps) {
  TestAdditiveGroup(zero, a, b, c, max_ulps);
  TestNonAbelianMultiplicativeGroup(one, c, x, y, max_ulps);
  TestVectorSpace(zero, a, b, c, zero, one, x, y, max_ulps);
}

}  // namespace testing_utilities
}  // namespace principia
