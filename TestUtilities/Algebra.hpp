#pragma once

#include "TestUtilities.hpp"

namespace Principia {
namespace TestUtilities {

template<typename T>
void TestEquality(T const& low, T const& high) {
  LogLine(L"Testing equality on " + ToString(low) + L" ≠ " +
          ToString(high) + L"...");
  AssertTrue(low == low, L"low == low was false.");
  AssertTrue(high == high, L"high == high was false.");
  AssertTrue(high != low, L"high != low was false.");
  AssertTrue(low != high, L"low != high was false.");

  LogLine(L"> True comparisons passed!");

  AssertTrue(!(high == low), L"high == low was true.");
  AssertTrue(!(low == high), L"low == high was true.");
  AssertTrue(!(low != low), L"low != low was true.");
  AssertTrue(!(high != high), L"high != high was true.");

  LogLine(L"> False comparisons passed!");
}

template<typename T>
void TestOrder(T const& low, T const& high) {
  TestEquality(low, high);

  LogLine(L"Testing ordering of " + ToString(low) + L" < " +
          ToString(high) + L"...");
  AssertTrue(high > low, L"high > low was false.");
  AssertTrue(low < high, L"low < high was false.");
  AssertTrue(low >= low, L"low >= low was false.");
  AssertTrue(low <= low, L"low <= low was false.");
  AssertTrue(high >= high, L"high >= high was false.");
  AssertTrue(high <= high, L"high <= high was false.");
  AssertTrue(high >= low, L"high >= low was false.");
  AssertTrue(low <= high, L"low <= high was false.");

  LogLine(L"> True comparisons passed!");

  AssertTrue(!(low > low), L"low > low was true.");
  AssertTrue(!(low < low), L"low < low was true.");
  AssertTrue(!(high > high), L"high > high was true.");
  AssertTrue(!(high < high), L"high < high was true.");
  AssertTrue(!(low > high), L"low > high was true.");
  AssertTrue(!(high < low), L"high < low was true.");
  AssertTrue(!(low >= high), L"low >= high was true.");
  AssertTrue(!(high <= low), L"high <= low was true.");

  LogLine(L"> False comparisons passed!");
}

template<typename T>
void TestAdditiveGroup(T const& zero, T const& a, T const& b, T const& c) {
  AssertEqual(a + zero, a);
  AssertEqual(zero + b, b);
  AssertEqual(a - a, zero);
  AssertEqual(-a - b, -(a + b));
  AssertEqual((a + b) + c, a + (b + c));
  AssertEqual(a - b - c, a - (b + c));
  AssertEqual(a + b, b + a);
  T accumulator = zero;
  accumulator += a;
  accumulator += b;
  accumulator -= c;
  AssertEqual(accumulator, a + b - c);
}

template<typename T>
void TestMultiplicativeGroup(T const& one, T const& a, T const& b, T const& c) {
  AssertEqual(a * one, a);
  AssertEqual(one * b, b);
  AssertEqual(a / a, one);
  AssertEqual((1 / a) / b, 1 / (a * b));
  AssertEqual((a * b) * c, a * (b * c));
  AssertEqual(a / b / c, a / (b * c));
  AssertEqual(a * b, b * a);
  T accumulator = one;
  accumulator *= a;
  accumulator *= b;
  accumulator /= c;
  AssertEqual(accumulator, a * b / c);
}

template<typename Map, typename Scalar, typename U, typename V>
void TestBilinearMap(Map const& map, U const& u1, U const& u2, V const& v1,
                     V const& v2, Scalar const& λ) {
  AssertEqual(map(u1 + u2, v1), map(u1, v1) + map(u2, v1));
  AssertEqual(map(u1, v1 + v2), map(u1, v1) + map(u1, v2));
  AssertEqual(map(λ * u1, v1), map(u1, λ * v1));
  AssertEqual(λ * map(u1, v1), map(u1, λ * v1));
  AssertEqual(map(u2 * λ, v2), map(u2, v2 * λ));
  AssertEqual(map(u2, v2) * λ, map(u2 * λ, v2));
}

template<typename Map, typename Scalar, typename U>
void TestSymmetricBilinearMap(Map const& map, U const& u1, U const& u2,
                               U const& v1, U const& v2, Scalar const& λ) {
  TestBilinearMap(map, u1, u2, v1, v2, λ);
  AssertEqual(map(u1, v1), map(v1, u1));
  AssertEqual(map(u2, v2), map(v2, u2));
}

template<typename Map, typename Scalar, typename U>
void TestAlternatingBilinearMap(Map const& map, U const& u1, U const& u2,
                                U const& v1, U const& v2, Scalar const& λ) {
  TestBilinearMap(map, u1, u2, v1, v2, λ);
  auto zero = map(u1, u1) - map(u1, u1);
  AssertEqual(map(u1, u1), zero);
  AssertEqual(map(u2, u2), zero);
  AssertEqual(map(v1, v2), -map(v2, v1));
}

template<typename Vector, typename Scalar>
void TestVectorSpace(Vector const& nullVector, Vector const& u, Vector const& v,
                     Vector const& w, Scalar const& zero, Scalar const& unit,
                     Scalar const& α, Scalar const& β) {
  TestAdditiveGroup(nullVector, u, v, w);
  AssertEqual((α * β) * v, α * (β * v));
  AssertEqual(unit * w, w);
  AssertEqual(u / unit, u);
  AssertEqual(zero * u, nullVector);
  AssertEqual(β * (u + v), β * u + v * β);
  AssertEqual((w + v) / α, w / α + v / α);
  AssertEqual((α + β) * w, α * w + β * w);
  AssertEqual(v * (α + β), α * v + β * v);
  Vector vector = u;
  vector *= α;
  AssertEqual(α * u, vector);
  vector /= α;
  AssertEqual(u, vector);
  vector *= zero;
  AssertEqual(vector, nullVector);
}

template<typename Vector, typename Scalar, typename Map>
void TestInnerProductSpace(Map const& map, Vector const& nullVector,
                           Vector const& u, Vector const& v, Vector const& w,
                           Vector const& a, Scalar const& zero,
                           Scalar const& unit, Scalar const& α,
                           Scalar const& β) {
  TestVectorSpace(nullVector, u, v, w, zero, unit, α, β);
  TestSymmetricBilinearMap(map, u, v, w, a, α);
}

template<typename T>
void TestField(T const& zero, T const& one, T const& a, T const& b,
               T const& c, T const& x, T const& y) {
  TestAdditiveGroup(zero, a, b, c);
  TestMultiplicativeGroup(one, c, x, y);
  TestVectorSpace(zero, a, b, c, zero, one, x, y);
}

}  // namespace TestUtilities
}  // namespace Principia
