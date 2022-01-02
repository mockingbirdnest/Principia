
#include "numerics/matrix_computations.hpp"

#include <tuple>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
    ::testing::Types<std::tuple<UnboundedVector<double>,
                                UnboundedLowerTriangularMatrix<double>,
                                UnboundedUpperTriangularMatrix<double>>>;

TYPED_TEST_SUITE(MatrixComputationsTest, MatrixTypes);

TYPED_TEST(MatrixComputationsTest, CholeskyDecomposition) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
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

  UpperTriangularMatrix r4_actual(4);
  Vector d4_actual(4);
  ᵗRDRDecomposition(hilbert4, r4_actual, d4_actual);
  EXPECT_THAT(d4_actual, AlmostEquals(d4_expected, 1615));
  EXPECT_THAT(r4_actual, AlmostEquals(r4_expected, 23));
}

TYPED_TEST(MatrixComputationsTest, BackSubstitution) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using UpperTriangularMatrix = typename std::tuple_element<2, TypeParam>::type;

  UpperTriangularMatrix const m3({1, 3, -2,
                                     4,  7,
                                         5});
  Vector const b3({1, 1, -4});
  Vector const x3_expected({-111.0 / 20.0, 33.0 / 20.0, -4.0 / 5.0});

  auto const x3_actual = BackSubstitution(m3, b3);
  EXPECT_THAT(x3_actual, AlmostEquals(x3_expected, 1));
}

TYPED_TEST(MatrixComputationsTest, ForwardSubstitution) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using LowerTriangularMatrix = typename std::tuple_element<1, TypeParam>::type;

  LowerTriangularMatrix const m3({1,
                                  3, -2,
                                  4,  7, 5});
  Vector const b3({1, 1, -4});
  Vector const x3_expected({1, 1, -3});

  auto const x3_actual = ForwardSubstitution(m3, b3);
  EXPECT_THAT(x3_actual, AlmostEquals(x3_expected, 0));
}

}  // namespace numerics
}  // namespace principia
