#include "numerics/matrix_views.hpp"

#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/unbounded_arrays.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_unbounded_arrays;

class MatrixViewsTest : public ::testing::Test {
 protected:
  MatrixViewsTest()
    : fu3_({6, -1, 12}),
      fv3_({10, 31, -47}),
      fv4_({-3, -3, 1, 4}),
      fm34_({-8,  -6, -4, -7,
             -4, -10,  9, -5,
              6,  -3, -2, -9}),
      fm23_({1, -2,  0,
             2,  3,  7}),
      fn23_({ 5, -1,  3,
             12, 13, -4}),
      uu3_({6, -1, 12}),
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

  FixedVector<double, 3> fu3_;
  FixedVector<double, 3> fv3_;
  FixedVector<double, 4> fv4_;
  FixedMatrix<double, 3, 4> fm34_;
  FixedMatrix<double, 2, 3> fm23_;
  FixedMatrix<double, 2, 3> fn23_;
  UnboundedVector<double> uu3_;
  UnboundedVector<double> uv3_;
  UnboundedVector<double> uv4_;
  UnboundedMatrix<double> um34_;
  UnboundedMatrix<double> um23_;
  UnboundedMatrix<double> un23_;
};

TEST_F(MatrixViewsTest, Indexing) {
  BlockView<FixedMatrix<double, 3, 4>> bfm34{.matrix = fm34_,
                                             .first_row = 1,
                                             .last_row = 2,
                                             .first_column = 0,
                                             .last_column = 2};
  EXPECT_EQ(-4, bfm34(0, 0));
  EXPECT_EQ(-10, bfm34(0, 1));
  EXPECT_EQ(6, bfm34(1, 0));
  EXPECT_EQ(-3, bfm34(1, 1));
}

}  // namespace numerics
}  // namespace principia
