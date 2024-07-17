#include "numerics/matrix_computations.hpp"

#include <tuple>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/transposed_view.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;
using namespace principia::testing_utilities::_vanishes_before;

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

TYPED_TEST(MatrixComputationsTest, ClassicalGramSchmidt) {
  using Matrix = typename std::tuple_element<3, TypeParam>::type;
  Matrix const m4({1, 2, 3, -4,
                   5, 6, 7, 8,
                   9, 8, -7, 6,
                   5, 4, 3, 2});
  auto const qr = ClassicalGramSchmidt(m4);

  // Check that the decomposition is correct.
  auto const near_m4 = qr.Q * Matrix(qr.R);
  EXPECT_THAT(near_m4, AlmostEquals(m4, 1));

  // Check that Q is nearly orthogonal.
  auto const near_identity = Matrix(TransposedView{.transpose = qr.Q}) * qr.Q;
  for (int i = 0; i < near_identity.rows(); ++i) {
    for (int j = 0; j < near_identity.columns(); ++j) {
      if (i == j) {
        EXPECT_THAT(near_identity(i, j), AlmostEquals(1, 0, 2));
      } else {
        EXPECT_THAT(near_identity(i, j), VanishesBefore(1, 0, 3));
      }
    }
  }
}

TYPED_TEST(MatrixComputationsTest, UnitriangularGramSchmidt) {
  using Matrix = typename std::tuple_element<3, TypeParam>::type;
  Matrix const m4({1, 2, 3, -4,
                   5, 6, 7, 8,
                   9, 8, -7, 6,
                   5, 4, 3, 2});
  auto const qr = UnitriangularGramSchmidt(m4);

  // Check that the decomposition is correct.
  auto const near_m4 = qr.Q * Matrix(qr.R);
  EXPECT_THAT(near_m4, AlmostEquals(m4, 1));

  // Check that R is unitriangular.
  for (int i = 0; i < qr.R.rows(); ++i) {
    EXPECT_THAT(qr.R(i, i), AlmostEquals(1, 0));
  }

  // Check that the columns of Q are approximately orthogonal.
  for (int c1 = 0; c1 < qr.Q.columns(); ++c1) {
    auto const column_c1 = ColumnView{.matrix = qr.Q,
                                      .first_row = 0,
                                      .last_row = qr.Q.rows() - 1,
                                      .column = c1};
    for (int c2 = c1 + 1; c2 < qr.Q.columns(); ++c2) {
      auto const column_c2 = ColumnView{.matrix = qr.Q,
                                        .first_row = 0,
                                        .last_row = qr.Q.rows() - 1,
                                        .column = c2};
      EXPECT_THAT(TransposedView{.transpose = column_c1} * column_c2,  // NOLINT
                  VanishesBefore(1, 24, 176));
    }
  }
}

TYPED_TEST(MatrixComputationsTest, UnitriangularGramSchmidt_Singular) {
  using Matrix = typename std::tuple_element<3, TypeParam>::type;
  Matrix const m4({1, 2, 0, -4,
                   5, 6, 0, 8,
                   9, 8, 0, 6,
                   5, 4, 0, 2});
  auto const qr = UnitriangularGramSchmidt(m4);

  // Check that the decomposition is correct.
  auto const near_m4 = qr.Q * Matrix(qr.R);
  EXPECT_THAT(near_m4, AlmostEquals(m4, 0));

  // Check that R is unitriangular.
  for (int i = 0; i < qr.R.rows(); ++i) {
    EXPECT_THAT(qr.R(i, i), AlmostEquals(1, 0));
  }

  // Check that the columns of Q are approximately orthogonal.
  for (int c1 = 0; c1 < qr.Q.columns(); ++c1) {
    auto const column_c1 = ColumnView{.matrix = qr.Q,
                                      .first_row = 0,
                                      .last_row = qr.Q.rows() - 1,
                                      .column = c1};
    for (int c2 = c1 + 1; c2 < qr.Q.columns(); ++c2) {
      auto const column_c2 = ColumnView{.matrix = qr.Q,
                                        .first_row = 0,
                                        .last_row = qr.Q.rows() - 1,
                                        .column = c2};
      EXPECT_THAT(TransposedView{.transpose = column_c1} * column_c2,  // NOLINT
                  VanishesBefore(1, 0, 64));
    }
  }
}

TYPED_TEST(MatrixComputationsTest, HessenbergForm) {
  using Matrix = typename std::tuple_element<3, TypeParam>::type;
  Matrix const m4({1, 2, 3, -4,
                   5, 6, 7, 8,
                   9, 8, -7, 6,
                   5, 4, 3, 2});
  Matrix const h4_expected({  1,
                              1.4852968963237645012,
                             -0.63177627362517020518,
                              5.1375822981102002423,

                             11.445523142259597039,
                              7.7328244274809160305,
                             10.840153034634636509,
                             -4.4908563773743662120,

                              0,
                             10.020455740333542876,
                             -3.9158079318984636532,
                              0.14244486136993501485,

                              0,
                              0,
                              2.4140754086886336638,
                             -2.8170164955824523773});
  auto h4_actual = HessenbergDecomposition(m4).H;
  // This component should really use VanishesBefore, but we don't have a good
  // way to do that.
  EXPECT_THAT(h4_actual(3, 1), VanishesBefore(1, 24));
  h4_actual(3, 1) = 0;
  EXPECT_THAT(h4_actual, AlmostEquals(h4_expected, 14));
}

TYPED_TEST(MatrixComputationsTest, RealSchurDecomposition) {
  using Matrix = typename std::tuple_element<3, TypeParam>::type;
  Matrix const m4({ 5,  4, -1,  0,
                    8, -1,  9,  8,
                   -4, -7,  2, -7,
                    8, -9, -2,  4});
  auto s4 = RealSchurDecomposition(m4, 1e-6);
  // Only check the real eigenvalues.
  EXPECT_THAT(
      s4.real_eigenvalues,
      ElementsAre(
          RelativeErrorFrom(6.2103405225078473234, IsNear(8.4e-7_(1))),
          RelativeErrorFrom(8.8004352424313246181, IsNear(6.0e-7_(1)))));
}

TYPED_TEST(MatrixComputationsTest, ClassicalJacobi) {
  using Vector = typename std::tuple_element<0, TypeParam>::type;
  using Matrix = typename std::tuple_element<3, TypeParam>::type;

  Matrix const m4({ 1, 0, -2, 3,
                    0, 4,  8, 1,
                   -2, 8,  3, 5,
                    3, 1,  5, 2});

  auto const actual = ClassicalJacobi(m4, /*max_iterations=*/20);
  EXPECT_THAT(actual.eigenvalues,
              AlmostEquals(Vector({ 4.1113216733296883337,
                                   13.142723660386208470,
                                   -6.7936661606326601037,
                                   -0.46037917308323670016}),
                           5));
  EXPECT_THAT(actual.rotation,
              AlmostEquals(Matrix({ 0.69703188150800443363259536246,
                                   -0.0242745668874571481045489833007,
                                    0.348414169206928869401042467771,
                                   -0.62623068294334091553061725607,

                                   -0.253480622026314290230470255713,
                                    0.63670280841355494645059139667,
                                   -0.461338172568875015343083925357,
                                   -0.56349285580764520595180494534,

                                   -0.086667079266903233811093219163,
                                    0.68301357119387773794045748748,
                                    0.67900055330874283257183529919,
                                    0.254832351836950053168269666832,

                                    0.66511874713460630976547141298,
                                    0.357089261565933987671587390098,
                                   -0.452474204605165882087767104452,
                                    0.474732983530015631333696554458}),
                           4));
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
              AlmostEquals(1.61073818969894450774101079322, 1, 2));
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
