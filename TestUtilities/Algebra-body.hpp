#pragma once

#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace test_utilities {

template<typename T>
void TestEquality(T const& low, T const& high) {
  LogLine("Testing equality on " + ToString(low) + " != " +
          ToString(high) + "...");
  AssertTrue(low == low, "low == low was false.");
  AssertTrue(high == high, "high == high was false.");
  AssertTrue(high != low, "high != low was false.");
  AssertTrue(low != high, "low != high was false.");

  AssertFalse(high == low, "high == low was true.");
  AssertFalse(low == high, "low == high was true.");
  AssertFalse(low != low, "low != low was true.");
  AssertFalse(high != high, "high != high was true.");
}

template<typename T>
void TestOrder(T const& low, T const& high) {
  TestEquality(low, high);
  LogLine("Testing ordering of " + ToString(low) + " < " +
          ToString(high) + "...");
  AssertTrue(high > low, "high > low was false.");
  AssertTrue(low < high, "low < high was false.");
  AssertTrue(low >= low, "low >= low was false.");
  AssertTrue(low <= low, "low <= low was false.");
  AssertTrue(high >= high, "high >= high was false.");
  AssertTrue(high <= high, "high <= high was false.");
  AssertTrue(high >= low, "high >= low was false.");
  AssertTrue(low <= high, "low <= high was false.");

  AssertFalse(low > low, "low > low was true.");
  AssertFalse(low < low, "low < low was true.");
  AssertFalse(high > high, "high > high was true.");
  AssertFalse(high < high, "high < high was true.");
  AssertFalse(low > high, "low > high was true.");
  AssertFalse(high < low, "high < low was true.");
  AssertFalse(low >= high, "low >= high was true.");
  AssertFalse(high <= low, "high <= low was true.");
}

template<typename T>
void TestAdditiveGroup(T const& zero, T const& a, T const& b, T const& c,
                       quantities::Dimensionless const& ε) {
  AssertEqual(a, +a);
  AssertEqual(a + zero, a);
  AssertEqual(zero + b, b);
  AssertEqual(a - a, zero);
  AssertEqual(-a - b, -(a + b));
  AssertEqual((a + b) + c, a + (b + c), ε);
  AssertEqual(a - b - c, a - (b + c), ε);
  AssertEqual(a + b, b + a);
  T accumulator = zero;
  accumulator += a;
  accumulator += b;
  accumulator -= c;
  AssertEqual(accumulator, a + b - c, ε);
}

template<typename T>
void TestMultiplicativeGroup(T const& one, T const& a, T const& b, T const& c,
                             quantities::Dimensionless const& ε) {
  AssertEqual(a * one, a);
  AssertEqual(one * b, b);
  AssertEqual(a / a, one);
  AssertEqual((1 / a) / b, 1 / (a * b), ε);
  AssertEqual((a * b) * c, a * (b * c), ε);
  AssertEqual(a / b / c, a / (b * c), ε);
  AssertEqual(a * b, b * a);
  T accumulator = one;
  accumulator *= a;
  accumulator *= b;
  accumulator /= c;
  AssertEqual(accumulator, a * b / c, ε);
}

template<typename Map, typename Scalar, typename U, typename V>
void TestBilinearMap(Map const& map, U const& u1, U const& u2, V const& v1,
                     V const& v2, Scalar const& λ,
                     quantities::Dimensionless const& ε) {
  AssertEqual(map(u1 + u2, v1), map(u1, v1) + map(u2, v1), ε);
  AssertEqual(map(u1, v1 + v2), map(u1, v1) + map(u1, v2), ε);
  AssertEqual(map(λ * u1, v1), map(u1, λ * v1), ε);
  AssertEqual(λ * map(u1, v1), map(u1, λ * v1), ε);
  AssertEqual(map(u2 * λ, v2), map(u2, v2 * λ), ε);
  AssertEqual(map(u2, v2) * λ, map(u2 * λ, v2), ε);
}

template<typename Map, typename Scalar, typename U>
void TestSymmetricBilinearMap(Map const& map, U const& u1, U const& u2,
                               U const& v1, U const& v2, Scalar const& λ,
                               quantities::Dimensionless const& ε) {
  TestBilinearMap(map, u1, u2, v1, v2, λ, ε);
  AssertEqual(map(u1, v1), map(v1, u1), ε);
  AssertEqual(map(u2, v2), map(v2, u2), ε);
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
  AssertTrue(map(u1, u1) > zero);
  AssertTrue(map(u2, u2) > zero);
  AssertTrue(map(v1, v1) > zero);
  AssertTrue(map(v2, v2) > zero);
}

template<typename Map, typename Scalar, typename U>
void TestAlternatingBilinearMap(Map const& map, U const& u1, U const& u2,
                                U const& v1, U const& v2, Scalar const& λ,
                                quantities::Dimensionless const& ε) {
  TestBilinearMap(map, u1, u2, v1, v2, λ, ε);
  auto zero = map(u1, u1) - map(u1, u1);
  AssertEqual(map(u1, u1), zero, ε);
  AssertEqual(map(u2, u2), zero, ε);
  AssertEqual(map(v1, v2), -map(v2, v1), ε);
}

template<typename Map, typename Scalar, typename U>
void TestLieBracket(Map const& map, U const& u1, U const& u2, U const& v1,
                    U const& v2, Scalar const& λ,
                    quantities::Dimensionless const& ε) {
  TestAlternatingBilinearMap(map, u1, u2, v1, v2, λ, ε);
  auto zero = map(u1, u1) - map(u1, u1);
  AssertEqual(map(u1, map(u2, v1)) +
              map(u2, map(v1, u1)) +
              map(v1, map(u1, u2)), zero, ε);
}

template<typename Vector, typename Scalar>
void TestVectorSpace(Vector const& nullVector, Vector const& u, Vector const& v,
                     Vector const& w, Scalar const& zero, Scalar const& unit,
                     Scalar const& α, Scalar const& β,
                     quantities::Dimensionless const& ε) {
  TestAdditiveGroup(nullVector, u, v, w, ε);
  AssertEqual((α * β) * v, α * (β * v), ε);
  AssertEqual(unit * w, w);
  AssertEqual(u / unit, u);
  AssertEqual(zero * u, nullVector);
  AssertEqual(β * (u + v), β * u + v * β, ε);
  AssertEqual((w + v) / α, w / α + v / α, ε);
  AssertEqual((α + β) * w, α * w + β * w, ε);
  AssertEqual(v * (α + β), α * v + β * v, ε);
  Vector vector = u;
  vector *= α;
  AssertEqual(α * u, vector);
  vector /= α;
  AssertEqual(u, vector, ε);
  vector *= zero;
  AssertEqual(vector, nullVector);
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

}  // namespace test_utilities
}  // namespace principia
