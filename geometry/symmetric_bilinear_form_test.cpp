
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
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Length;
using quantities::Pow;
using quantities::Square;
using quantities::si::Metre;
using testing_utilities::EqualsProto;
using ::testing::Eq;

class SymmetricBilinearFormTest : public ::testing::Test {
 protected:
  using World =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST, true>;
  using Eigenworld =
    Frame<serialization::Frame::TestTag, serialization::Frame::TEST1, false>;

  SymmetricBilinearFormTest() {}

  SymmetricBilinearForm<Length, World> MakeSymmetricBilinearForm(
      R3x3Matrix<double> const& untyped_matrix) {
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

TEST_F(SymmetricBilinearFormTest, Equality) {
  auto const f1 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                               {2, -3,  5},
                                                               {4,  5,  0}));
  auto const f2 = MakeSymmetricBilinearForm(R3x3Matrix<double>({1, 2, 3},
                                                               {2, 4, 0},
                                                               {3, 0, 5}));
  EXPECT_TRUE(f1 == f1);
  EXPECT_TRUE(f1 != f2);
}

TEST_F(SymmetricBilinearFormTest, UnaryOperators) {
  auto const f = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}));
  EXPECT_THAT(+f,
              Eq(MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}))));
  EXPECT_THAT(-f,
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

TEST_F(SymmetricBilinearFormTest, LinearMap) {
  auto const f = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}));
  Vector<Length, World> v({1.0 * Metre, 3.0 * Metre, -1.0 * Metre});
  EXPECT_THAT(
      f * v,
      Eq(Vector<Square<Length>, World>({3.0 * Pow<2>(Metre),
                                        -12.0 * Pow<2>(Metre),
                                        19.0 * Pow<2>(Metre)})));
  EXPECT_THAT(
      v * f,
      Eq(Vector<Square<Length>, World>({3.0 * Pow<2>(Metre),
                                        -12.0 * Pow<2>(Metre),
                                        19.0 * Pow<2>(Metre)})));
}

TEST_F(SymmetricBilinearFormTest, SymmetricProduct) {
  Vector<Length, World> v1({1.0 * Metre, 3.0 * Metre, -1.0 * Metre});
  Vector<double, World> v2({2.0, 6.0, -5.0});
  EXPECT_THAT(SymmetricProduct(v1, v2),
              Eq(MakeSymmetricBilinearForm(
                     R3x3Matrix<double>({   2,     6,  -3.5},
                                        {   6,    18, -10.5},
                                        {-3.5, -10.5,   5}))));
}

TEST_F(SymmetricBilinearFormTest, Anticommutator) {
  auto const f = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}));
  Bivector<Length, World> const b({1.0 * Metre, 3.0 * Metre, -5.0 * Metre});
  EXPECT_THAT(
      Anticommutator(f, b),
      Eq(Bivector<Square<Length>, World>({11 * Pow<2>(Metre),
                                          26 * Pow<2>(Metre),
                                          -9 * Pow<2>(Metre)})));
}

TEST_F(SymmetricBilinearFormTest, InnerProductForm) {
  Vector<Length, World> v1({1.0 * Metre, 3.0 * Metre, -1.0 * Metre});
  Vector<double, World> v2({2.0, 6.0, -5.0});
  auto const a = InnerProductForm<World>()(v1, v2);
  EXPECT_THAT(InnerProductForm<World>()(v1, v2), Eq(25 * Metre));
}

TEST_F(SymmetricBilinearFormTest, Apply) {
  auto const f = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}));
  Vector<Length, World> v1({1.0 * Metre, 3.0 * Metre, -1.0 * Metre});
  Vector<Length, World> v2({2.0 * Metre, 6.0 * Metre, -5.0 * Metre});
  EXPECT_THAT(f(v1, v2), Eq(-161 * Pow<3>(Metre)));
}

TEST_F(SymmetricBilinearFormTest, Serialization) {
  auto const f = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  2,  4},
                                                              {2, -3,  5},
                                                              {4,  5,  0}));
  serialization::SymmetricBilinearForm message1;
  f.WriteToMessage(&message1);
  EXPECT_TRUE(message1.has_frame());
  EXPECT_TRUE(message1.has_matrix());
  auto const g =
      SymmetricBilinearForm<Length, World>::ReadFromMessage(message1);
  EXPECT_EQ(f, g);
  serialization::SymmetricBilinearForm message2;
  g.WriteToMessage(&message2);
  EXPECT_THAT(message2, EqualsProto(message1));
}

TEST_F(SymmetricBilinearFormTest, Diagonalize) {
  auto const f = MakeSymmetricBilinearForm(R3x3Matrix<double>({1,  0,  0},
                                                              {0, -3,  0},
                                                              {0,  0,  2}));
  auto const fd = f.Diagonalize<Eigenworld>();
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
