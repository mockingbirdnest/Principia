
#include "numerics/matrix_computations.hpp"

#include <tuple>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using quantities::Sqrt;
using testing_utilities::AlmostEquals;

template<typename T>
class MatrixComputationsTest : public ::testing::Test {
 protected:
};

using MatrixTypes =
    ::testing::Types<std::tuple<FixedVector<double, 4>,
                                FixedLowerTriangularMatrix<double, 4>,
                                FixedUpperTriangularMatrix<double, 4>,
                                FixedMatrix<double, 4, 4>>,
                     std::tuple<UnboundedVector<double>,
                                UnboundedLowerTriangularMatrix<double>,
                                UnboundedUpperTriangularMatrix<double>,
                                UnboundedMatrix<double>>>;

TYPED_TEST_SUITE(MatrixComputationsTest, MatrixTypes);

TYPED_TEST(MatrixComputationsTest, CholeskyDecomposition) {
  using UpperTriangularMatrix = typename std::tuple_element<2, TypeParam>::type;

  UpperTriangularMatrix const hilbert4({1, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0,
                                           1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0,
                                                      1.0 / 5.0, 1.0 / 6.0,
                                                                 1.0 / 7.0});
  UpperTriangularMatrix const r4_expected({
      1,        1.0 / 2.0,         1.0 / 3.0,          1.0 / 4.0,
         1.0 / Sqrt(12.0),  1.0 / Sqrt(12.0),  Sqrt(27.0) / 20.0,
                           1.0 / Sqrt(180.0),   1.0 / Sqrt(80.0),
                                              1.0 / Sqrt(2800.0)});

  auto const r4_actual = CholeskyDecomposition(hilbert4);
  EXPECT_THAT(r4_actual, AlmostEquals(r4_expected, 245));
}

TYPED_TEST(MatrixComputationsTest, ᵗRDRDecomposition) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using UpperTriangularMatrix = typename std::tuple_element<2, TypeParam>::type;

  UpperTriangularMatrix const hilbert4({1, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0,
                                           1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0,
                                                      1.0 / 5.0, 1.0 / 6.0,
                                                                 1.0 / 7.0});
  UpperTriangularMatrix const r4_expected({
      1, 1.0 / 2.0, 1.0 / 3.0,  1.0 / 4.0,
                 1,         1, 9.0 / 10.0,
                            1,  3.0 / 2.0,
                                        1});
  Vector d4_expected({1, 1.0 / 12.0, 1.0 / 180.0, 1.0 / 2800.0});

  const auto [r4_actual, d4_actual] = ᵗRDRDecomposition<Vector>(hilbert4);
  EXPECT_THAT(d4_actual, AlmostEquals(d4_expected, 1615));
  EXPECT_THAT(r4_actual, AlmostEquals(r4_expected, 23));
}

TYPED_TEST(MatrixComputationsTest, BackSubstitution) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using UpperTriangularMatrix = typename std::tuple_element<2, TypeParam>::type;

  UpperTriangularMatrix const m4({1, 3, -2,  6,
                                     4,  7, -1,
                                         5,  3,
                                             2});
  Vector const b4({1, 1, -4, 4});
  Vector const x4_expected({-111.0 / 4.0, 17.0 / 4.0, -2.0, 2.0});

  auto const x4_actual = BackSubstitution(m4, b4);
  EXPECT_THAT(x4_actual, AlmostEquals(x4_expected, 0));
}

TYPED_TEST(MatrixComputationsTest, ForwardSubstitution) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using LowerTriangularMatrix = typename std::tuple_element<1, TypeParam>::type;

  LowerTriangularMatrix const m4({ 1,
                                   3, -2,
                                   4,  7, 5,
                                  -6,  9, 1, 2});
  Vector const b4({1, 1, -4, 4});
  Vector const x4_expected({1, 1, -3, 2});

  auto const x4_actual = ForwardSubstitution(m4, b4);
  EXPECT_THAT(x4_actual, AlmostEquals(x4_expected, 0));
}

TYPED_TEST(MatrixComputationsTest, RayleighQuotient) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using Matrix = typename std::tuple_element<3, TypeParam>::type;

  Matrix const m4({ 1, 0, -2, 3,
                   -4, 4,  1, 2,
                    0, 8,  3, 5,
                   -7, 1,  2, 2});
  Vector const v4({1, -1, 2, 3});

  auto const actual = RayleighQuotient(m4, v4);
  EXPECT_THAT(actual, AlmostEquals(38.0 / 15.0, 0));
}

TYPED_TEST(MatrixComputationsTest, RayleighQuotientIteration) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using Matrix = typename std::tuple_element<3, TypeParam>::type;

  Matrix const m4({ 1, 0, -2, 3,
                   -4, 4,  1, 2,
                    0, 8,  3, 5,
                   -7, 1,  2, 2});
  Vector const v4({1, -1, 2, 3});

  auto const actual = RayleighQuotientIteration(m4, v4);
  EXPECT_THAT(actual.eigenvalue,
              AlmostEquals(1.61073818969894450774101079322, 1));
  EXPECT_THAT(actual.eigenvector,
              AlmostEquals(Vector({-0.169506517685592297780621340523,
                                   0.444078911480267347768558852121,
                                   -0.71572747951425189504245569810,
                                   -0.51165968759399151699548201390}),
                           1));
}

TYPED_TEST(MatrixComputationsTest, Solve) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using Matrix = typename std::tuple_element<3, TypeParam>::type;

  Matrix const m4({ 1, 0, -2, 3,
                   -4, 4,  1, 2,
                    0, 8,  3, 5,
                   -7, 1,  2, 2});
  Vector const v4({1, -1, 2, 3});
  Vector const x4_expected(
      {66.0 / 383.0, -397.0 / 383.0, 539.0 / 383.0, 465.0 / 383.0});

  auto const x4_actual = Solve(m4, v4);
  EXPECT_THAT(x4_actual, AlmostEquals(x4_expected, 4));
}

}  // namespace numerics
}  // namespace principia
