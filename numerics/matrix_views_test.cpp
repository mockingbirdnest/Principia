#include "numerics/matrix_views.hpp"

#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/transposed_view.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::testing_utilities::_almost_equals;

class MatrixViewsTest : public ::testing::Test {
 protected:
  MatrixViewsTest()
    : fu2_({6, -1}),
      fv3_({10, 31, -47}),
      fv4_({-3, -3, 1, 4}),
      fm34_({-8,  -6, -4, -7,
             -4, -10,  9, -5,
              6,  -3, -2, -9}),
      fm23_({1, -2,  0,
             2,  3,  7}),
      fn23_({ 5, -1,  3,
             12, 13, -4}),
      uu2_({6, -1}),
      uv3_({10, 31, -47}),
      uv4_({-3, -3, 1, 4}),
      um34_(3, 4,
            {-8,  -6, -4, -7,
             -4, -10,  9, -5,
              6,  -3, -2, -9}),
      um23_(2, 3,
            {1, -2,  0,
             2,  3,  7}),
      un23_(2, 3,
            { 5, -1,  3,
             12, 13, -4}) {}

  FixedVector<double, 2> fu2_;
  FixedVector<double, 3> fv3_;
  FixedVector<double, 4> fv4_;
  FixedMatrix<double, 3, 4> fm34_;
  FixedMatrix<double, 2, 3> fm23_;
  FixedMatrix<double, 2, 3> fn23_;
  UnboundedVector<double> uu2_;
  UnboundedVector<double> uv3_;
  UnboundedVector<double> uv4_;
  UnboundedMatrix<double> um34_;
  UnboundedMatrix<double> um23_;
  UnboundedMatrix<double> un23_;
};

TEST_F(MatrixViewsTest, BlockView_Indexing) {
  BlockView<FixedMatrix<double, 3, 4>> bfm34{.matrix = fm34_,
                                             .first_row = 1,
                                             .last_row = 2,
                                             .first_column = 0,
                                             .last_column = 2};
  EXPECT_EQ(-4, bfm34(0, 0));
  EXPECT_EQ(-10, bfm34(0, 1));
  EXPECT_EQ(6, bfm34(1, 0));
  EXPECT_EQ(-3, bfm34(1, 1));

  BlockView<UnboundedMatrix<double>> bum34{.matrix = um34_,
                                           .first_row = 1,
                                           .last_row = 2,
                                           .first_column = 0,
                                           .last_column = 2};
  EXPECT_EQ(-4, bum34(0, 0));
  EXPECT_EQ(-10, bum34(0, 1));
  EXPECT_EQ(6, bum34(1, 0));
  EXPECT_EQ(-3, bum34(1, 1));
}

TEST_F(MatrixViewsTest, BlockView_Addition) {
  BlockView<FixedMatrix<double, 3, 4>> bfm34{.matrix = fm34_,
                                             .first_row = 1,
                                             .last_row = 2,
                                             .first_column = 0,
                                             .last_column = 2};
  bfm34 += um23_;
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     {-3, -12, 9,
                                       8,   0, 5})),
             bfm34);

  BlockView<UnboundedMatrix<double>> bum34{.matrix = um34_,
                                           .first_row = 1,
                                           .last_row = 2,
                                           .first_column = 0,
                                           .last_column = 2};
  bum34 -= fm23_;
  EXPECT_EQ((FixedMatrix<double, 2, 3>({-5, -8,  9,
                                         4, -6, -9})),
             bum34);
}

TEST_F(MatrixViewsTest, BlockView_Multiplication) {
  BlockView<FixedMatrix<double, 3, 4>> bfm34{.matrix = fm34_,
                                             .first_row = 1,
                                             .last_row = 2,
                                             .first_column = 0,
                                             .last_column = 2};
  bfm34 *= 3;
  EXPECT_EQ((UnboundedMatrix<double>(2, 3,
                                     {-12, -30, 27,
                                       18,  -9, -6})),
             bfm34);

  BlockView<UnboundedMatrix<double>> bum34{.matrix = um34_,
                                           .first_row = 1,
                                           .last_row = 2,
                                           .first_column = 0,
                                           .last_column = 2};
  bum34 /= 2;
  EXPECT_EQ((FixedMatrix<double, 2, 3>({-2,   -5,  4.5,
                                         3, -1.5,   -1})),
             bum34);
}

TEST_F(MatrixViewsTest, ColumnView_Indexing) {
  ColumnView<FixedMatrix<double, 3, 4>> cfm34{.matrix = fm34_,
                                              .first_row = 1,
                                              .last_row = 2,
                                              .column = 3};
  EXPECT_EQ(-5, cfm34[0]);
  EXPECT_EQ(-9, cfm34[1]);

  ColumnView<UnboundedMatrix<double>> cum34{.matrix = um34_,
                                            .first_row = 1,
                                            .last_row = 2,
                                            .column = 3};
  EXPECT_EQ(-5, cum34[0]);
  EXPECT_EQ(-9, cum34[1]);
}

TEST_F(MatrixViewsTest, ColumnView_Assignment) {
  ColumnView<FixedMatrix<double, 3, 4>> cfm34{.matrix = fm34_,
                                              .first_row = 1,
                                              .last_row = 2,
                                              .column = 2};
  ColumnView<UnboundedMatrix<double>> cum34{.matrix = um34_,
                                            .first_row = 1,
                                            .last_row = 2,
                                            .column = 3};
  cum34 = cfm34;
  EXPECT_EQ(9, cum34[0]);
  EXPECT_EQ(-2, cum34[1]);

  cfm34 = cfm34;
}

TEST_F(MatrixViewsTest, ColumnView_Addition) {
  ColumnView<FixedMatrix<double, 3, 4>> cfm34{.matrix = fm34_,
                                              .first_row = 1,
                                              .last_row = 2,
                                              .column = 3};
  cfm34 += uu2_;
  EXPECT_EQ((UnboundedVector<double>({1, -10})), cfm34);

  ColumnView<UnboundedMatrix<double>> cum34{.matrix = um34_,
                                            .first_row = 1,
                                            .last_row = 2,
                                            .column = 3};
  cum34 -= fu2_;
  EXPECT_EQ((FixedVector<double, 2>({-11, -8})), cum34);
}

TEST_F(MatrixViewsTest, ColumnView_Multiplication) {
  ColumnView<FixedMatrix<double, 3, 4>> cfm34{.matrix = fm34_,
                                              .first_row = 1,
                                              .last_row = 2,
                                              .column = 3};
  cfm34 *= 3;
  EXPECT_EQ((UnboundedVector<double>({-15, -27})), cfm34);

  ColumnView<UnboundedMatrix<double>> cum34{.matrix = um34_,
                                            .first_row = 1,
                                            .last_row = 2,
                                            .column = 3};
  cum34 /= 2;
  EXPECT_EQ((FixedVector<double, 2>({-2.5, -4.5})), cum34);
}

TEST_F(MatrixViewsTest, ColumnView_Norm) {
  ColumnView<FixedMatrix<double, 3, 4>> cfm34{.matrix = fm34_,
                                              .first_row = 1,
                                              .last_row = 2,
                                              .column = 3};
  EXPECT_EQ(106, cfm34.Norm²());
  EXPECT_THAT(cfm34.Norm(), AlmostEquals(Sqrt(106), 0));

  ColumnView<UnboundedMatrix<double>> cum34{.matrix = um34_,
                                            .first_row = 1,
                                            .last_row = 2,
                                            .column = 3};
  EXPECT_EQ(106, cum34.Norm²());
  EXPECT_THAT(cum34.Norm(), AlmostEquals(Sqrt(106), 0));
}

TEST_F(MatrixViewsTest, ColumnView_DotProduct) {
  EXPECT_EQ(
      -60,
      (TransposedView{.transpose = ColumnView{.matrix = fm34_,  // NOLINT
                                              .first_row = 0,
                                              .last_row = 2,
                                              .column = 1}} *
       ColumnView{.matrix = um34_,
                  .first_row = 0,
                  .last_row = 2,
                  .column = 2}));
}

}  // namespace numerics
}  // namespace principia
