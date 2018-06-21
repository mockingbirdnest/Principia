
#include "geometry/symmetric_bilinear_form.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Length;
using quantities::Square;
using quantities::si::Metre;
using ::testing::Eq;

class SymmetricBilinearFormTest : public ::testing::Test {
 protected:
  using World =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;

  SymmetricBilinearFormTest() {}

  SymmetricBilinearForm<Length, World> MakeSymmetricBilinearForm(
      R3x3Matrix<double> const& untyped_matrix) {
    CHECK_EQ(untyped_matrix(0, 1), untyped_matrix(1, 0));
    CHECK_EQ(untyped_matrix(0, 2), untyped_matrix(2, 0));
    CHECK_EQ(untyped_matrix(1, 2), untyped_matrix(2, 1));
    R3x3Matrix<Length> const typed_matrix(
        R3Element<double>(
            untyped_matrix(0, 0), untyped_matrix(0, 1), untyped_matrix(0, 2)) *
            Metre,
        R3Element<double>(
            untyped_matrix(1, 0), untyped_matrix(1, 1), untyped_matrix(1, 2)) *
            Metre,
        R3Element<double>(
            untyped_matrix(2, 0), untyped_matrix(2, 1), untyped_matrix(2, 2)) *
            Metre);
    return SymmetricBilinearForm<Length, World>(typed_matrix);
  }
};

TEST_F(SymmetricBilinearFormTest, UnaryOperators) {
  auto const f1 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                               {2, -3,  5},
                                                               {4,  5,  0}));
  EXPECT_THAT(+f1,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}))));
  EXPECT_THAT(-f1,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({-1, -2, -4},
                                                              {-2,  3, -5},
                                                              {-4, -5,  0}))));
}

TEST_F(SymmetricBilinearFormTest, BinaryOperators) {
  auto const f1 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                               {2, -3,  5},
                                                               {4,  5,  0}));
  auto const f2 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1, 2, 3},
                                                               {2, 4, 0},
                                                               {3, 0, 5}));
  EXPECT_THAT(f1 + f2,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({2, 4, 7},
                                                              {4, 1, 5},
                                                              {7, 5, 5}))));
  EXPECT_THAT(f1 - f2,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({0,  0,  1},
                                                              {0, -7,  5},
                                                              {1,  5,  -5}))));
  EXPECT_THAT((-2) * f1,
              Eq(MakeSymmetricBilinearForm(
                     R3x3Matrix<double>({-2,  -4,  -8},
                                        {-4,   6, -10},
                                        {-8, -10,   0}))));
  EXPECT_THAT(f1 * 3,
              Eq(MakeSymmetricBilinearForm(
                     R3x3Matrix<double>({ 3,  6, 12},
                                        { 6, -9, 15},
                                        {12, 15,  0}))));
  EXPECT_THAT(f1 / 2,
              Eq(MakeSymmetricBilinearForm(
                     R3x3Matrix<double>({0.5,    1,   2},
                                        {  1, -1.5, 2.5},
                                        {  2,  2.5,   0}))));
}

TEST_F(SymmetricBilinearFormTest, Assignment) {
  auto const f1 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                               {2, -3,  5},
                                                               {4,  5,  0}));
  auto const f2 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1, 2, 3},
                                                               {2, 4, 0},
                                                               {3, 0, 5}));
  auto a = f1;
  auto b = f1;
  auto c = f1;
  auto d = f1;
  a += f2;
  b -= f2;
  c *= 3;
  d /= 2;
  EXPECT_THAT(a,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({2, 4, 7},
                                                              {4, 1, 5},
                                                              {7, 5, 5}))));
  EXPECT_THAT(b,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({0,  0,  1},
                                                              {0, -7,  5},
                                                              {1,  5,  -5}))));
  EXPECT_THAT(c,
              Eq(MakeSymmetricBilinearForm(
                     R3x3Matrix<double>({ 3,  6, 12},
                                        { 6, -9, 15},
                                        {12, 15,  0}))));
  EXPECT_THAT(d,
              Eq(MakeSymmetricBilinearForm(
                     R3x3Matrix<double>({0.5,    1,   2},
                                        {  1, -1.5, 2.5},
                                        {  2,  2.5,   0}))));
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
