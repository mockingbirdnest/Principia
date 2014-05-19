#pragma once

#include "gtest/gtest.h"
#include "testing_utilities/algebra.hpp"

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
                       quantities::Dimensionless const& ε) {
  ASSERT_EQ(a, +a);
  ASSERT_EQ(a + zero, a);
  ASSERT_EQ(zero + b, b);
  ASSERT_EQ(a - a, zero);
  ASSERT_EQ(-a - b, -(a + b));
  ASSERT_EQ((a + b) + c, a + (b + c), ε);
  ASSERT_EQ(a - b - c, a - (b + c), ε);
  ASSERT_EQ(a + b, b + a);
  T accumulator = zero;
  accumulator += a;
  accumulator += b;
  accumulator -= c;
  ASSERT_EQ(accumulator, a + b - c, ε);
}

template<typename T>
void TestMultiplicativeGroup(T const& one, T const& a, T const& b, T const& c,
                             quantities::Dimensionless const& ε) {
  ASSERT_EQ(a * one, a);
  ASSERT_EQ(one * b, b);
  ASSERT_EQ(a / a, one);
  ASSERT_EQ((1 / a) / b, 1 / (a * b), ε);
  ASSERT_EQ((a * b) * c, a * (b * c), ε);
  ASSERT_EQ(a / b / c, a / (b * c), ε);
  ASSERT_EQ(a * b, b * a);
  T accumulator = one;
  accumulator *= a;
  accumulator *= b;
  accumulator /= c;
  ASSERT_EQ(accumulator, a * b / c, ε);
}

template<typename Map, typename Scalar, typename U, typename V>
void TestBilinearMap(Map const& map, U const& u1, U const& u2, V const& v1,
                     V const& v2, Scalar const& λ,
                     quantities::Dimensionless const& ε) {
  ASSERT_EQ(map(u1 + u2, v1), map(u1, v1) + map(u2, v1), ε);
  ASSERT_EQ(map(u1, v1 + v2), map(u1, v1) + map(u1, v2), ε);
  ASSERT_EQ(map(λ * u1, v1), map(u1, λ * v1), ε);
  ASSERT_EQ(λ * map(u1, v1), map(u1, λ * v1), ε);
  ASSERT_EQ(map(u2 * λ, v2), map(u2, v2 * λ), ε);
  ASSERT_EQ(map(u2, v2) * λ, map(u2 * λ, v2), ε);
}

template<typename Map, typename Scalar, typename U>
void TestSymmetricBilinearMap(Map const& map, U const& u1, U const& u2,
                               U const& v1, U const& v2, Scalar const& λ,
                               quantities::Dimensionless const& ε) {
  TestBilinearMap(map, u1, u2, v1, v2, λ, ε);
  ASSERT_EQ(map(u1, v1), map(v1, u1), ε);
  ASSERT_EQ(map(u2, v2), map(v2, u2), ε);
}

template<typename Map, typename Scalar, typename U>
void TestSymmetricPositiveDefiniteBilinearMap(
    Map const& map,
    U const& u1,
    U const& u2, U const& v1,
    U const& v2, Scalar const& λ,
    quantities::Dimensionless const& ε) {
  TestSymmetricBilinearMap(map, u1, u2, v1, v2, λ, ε);
  auto zero = map(u1, u1) - map(u1, u1);
  EXPECT_TRUE(map(u1, u1) > zero) << map(u1, u1);
  EXPECT_TRUE(map(u2, u2) > zero) << map(u2, u2);
  EXPECT_TRUE(map(v1, v1) > zero) << map(v1, v1);
  EXPECT_TRUE(map(v2, v2) > zero) << map(v2, v2);
}

template<typename Map, typename Scalar, typename U>
void TestAlternatingBilinearMap(Map const& map, U const& u1, U const& u2,
                                U const& v1, U const& v2, Scalar const& λ,
                                quantities::Dimensionless const& ε) {
  TestBilinearMap(map, u1, u2, v1, v2, λ, ε);
  auto zero = map(u1, u1) - map(u1, u1);
  ASSERT_EQ(map(u1, u1), zero, ε);
  ASSERT_EQ(map(u2, u2), zero, ε);
  ASSERT_EQ(map(v1, v2), -map(v2, v1), ε);
}

template<typename Map, typename Scalar, typename U>
void TestLieBracket(Map const& map, U const& u1, U const& u2, U const& v1,
                    U const& v2, Scalar const& λ,
                    quantities::Dimensionless const& ε) {
  TestAlternatingBilinearMap(map, u1, u2, v1, v2, λ, ε);
  auto zero = map(u1, u1) - map(u1, u1);
  ASSERT_EQ(map(u1, map(u2, v1)) + map(u2, map(v1, u1)) + map(v1, map(u1, u2)),
            zero, ε);
}

template<typename Vector, typename Scalar>
void TestVectorSpace(Vector const& nullVector, Vector const& u, Vector const& v,
                     Vector const& w, Scalar const& zero, Scalar const& unit,
                     Scalar const& α, Scalar const& β,
                     quantities::Dimensionless const& ε) {
  TestAdditiveGroup(nullVector, u, v, w, ε);
  ASSERT_EQ((α * β) * v, α * (β * v), ε);
  ASSERT_EQ(unit * w, w);
  ASSERT_EQ(u / unit, u);
  ASSERT_EQ(zero * u, nullVector);
  ASSERT_EQ(β * (u + v), β * u + v * β, ε);
  ASSERT_EQ((w + v) / α, w / α + v / α, ε);
  ASSERT_EQ((α + β) * w, α * w + β * w, ε);
  ASSERT_EQ(v * (α + β), α * v + β * v, ε);
  Vector vector = u;
  vector *= α;
  ASSERT_EQ(α * u, vector);
  vector /= α;
  ASSERT_EQ(u, vector, ε);
  vector *= zero;
  ASSERT_EQ(vector, nullVector);
}

template<typename Vector, typename Scalar, typename Map>
void TestInnerProductSpace(Map const& map, Vector const& nullVector,
                           Vector const& u, Vector const& v, Vector const& w,
                           Vector const& a, Scalar const& zero,
                           Scalar const& unit, Scalar const& α, Scalar const& β,
                           quantities::Dimensionless const& ε) {
  TestVectorSpace(nullVector, u, v, w, zero, unit, α, β, ε);
  TestSymmetricPositiveDefiniteBilinearMap(map, u, v, w, a, α, ε);
}

template<typename T>
void TestField(T const& zero, T const& one, T const& a, T const& b, T const& c,
               T const& x, T const& y, quantities::Dimensionless const& ε) {
  TestAdditiveGroup(zero, a, b, c, ε);
  TestMultiplicativeGroup(one, c, x, y, ε);
  TestVectorSpace(zero, a, b, c, zero, one, x, y, ε);
}

}  // namespace testing_utilities
}  // namespace principia
