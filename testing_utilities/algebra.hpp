#pragma once

#include <cstdint>

namespace principia {
namespace testing_utilities {
namespace _algebra {
namespace internal {

template<typename T>
void TestEquality(T const& low, T const& high);

template<typename T>
void TestOrder(T const& low, T const& high);

template<typename T>
void TestAdditiveGroup(T const& zero, T const& a, T const& b, T const& c,
                       std::int64_t min_ulps = 0,
                       std::int64_t max_ulps = 0);

template<typename T>
void TestGroup(T const& identity, T const& a, T const& b, T const& c,
               T (*operation)(T const&, T const&), T (*inverse)(T const&),
               std::int64_t min_ulps = 0,
               std::int64_t max_ulps = 0);

template<typename T>
void TestAbelianMultiplicativeGroup(
    T const& one, T const& a, T const& b, T const& c,
    std::int64_t max_ulps = 0);

template<typename T>
void TestNonAbelianMultiplicativeGroup(
    T const& one, T const& a, T const& b, T const& c,
    std::int64_t min_ulps = 0,
    std::int64_t max_ulps = 0);

template<typename Map, typename Scalar, typename U, typename V>
void TestBilinearMap(Map const& map, U const& u1, U const& u2, V const& v1,
                     V const& v2, Scalar const& λ,
                     std::int64_t min_ulps = 0,
                     std::int64_t max_ulps = 0);

template<typename Map, typename Scalar, typename U>
void TestSymmetricBilinearMap(Map const& map, U const& u1, U const& u2,
                              U const& v1, U const& v2, Scalar const& λ,
                              std::int64_t min_ulps = 0,
                              std::int64_t max_ulps = 0);

template<typename Map, typename Scalar, typename U>
void TestSymmetricPositiveDefiniteBilinearMap(
    Map const& map, U const& u1, U const& u2,
    U const& v1, U const& v2, Scalar const& λ,
    std::int64_t min_ulps = 0,
    std::int64_t max_ulps = 0);

template<typename Map, typename Scalar, typename U>
void TestAlternatingBilinearMap(Map const& map, U const& u1, U const& u2,
                                U const& v1, U const& v2, Scalar const& λ,
                                std::int64_t min_ulps = 0,
                                std::int64_t max_ulps = 0);

template<typename Map, typename Scalar, typename U>
void TestLieBracket(Map const& map, U const& u1, U const& u2,
                    U const& v1, U const& v2, Scalar const& λ,
                    std::int64_t min_ulps = 0,
                    std::int64_t max_ulps = 0);

template<typename Vector, typename Scalar>
void TestVectorSpace(Vector const& nullVector, Vector const& u, Vector const& v,
                     Vector const& w, Scalar const& zero, Scalar const& unit,
                     Scalar const& α, Scalar const& β,
                     std::int64_t min_ulps = 0,
                     std::int64_t max_ulps = 0);

template<typename Vector, typename Scalar, typename Map>
void TestInnerProductSpace(Map const& map, Vector const& nullVector,
                           Vector const& u, Vector const& v, Vector const& w,
                           Vector const& a, Scalar const& zero,
                           Scalar const& unit, Scalar const& α,
                           Scalar const& β,
                           std::int64_t min_ulps = 0,
                           std::int64_t max_ulps = 0);

template<typename T>
void TestField(T const& zero, T const& one,
               T const& a, T const& b, T const& c, T const& x, T const& y,
               std::int64_t min_ulps = 0,
               std::int64_t max_ulps = 0);

template<typename T>
void TestSkewField(T const& zero, T const& one,
                   T const& a, T const& b, T const& c, T const& x, T const& y,
                   std::int64_t min_ulps = 0,
                   std::int64_t max_ulps = 0);

}  // namespace internal

using internal::TestAbelianMultiplicativeGroup;
using internal::TestAdditiveGroup;
using internal::TestAlternatingBilinearMap;
using internal::TestBilinearMap;
using internal::TestEquality;
using internal::TestField;
using internal::TestGroup;
using internal::TestInnerProductSpace;
using internal::TestLieBracket;
using internal::TestNonAbelianMultiplicativeGroup;
using internal::TestOrder;
using internal::TestSkewField;
using internal::TestSymmetricBilinearMap;
using internal::TestSymmetricPositiveDefiniteBilinearMap;
using internal::TestVectorSpace;

}  // namespace _algebra
}  // namespace testing_utilities
}  // namespace principia

namespace principia::testing_utilities {
using namespace principia::testing_utilities::_algebra;
}  // namespace principia::testing_utilities

#include "testing_utilities/algebra_body.hpp"
