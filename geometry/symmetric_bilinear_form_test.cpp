
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
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Length;
using quantities::Pow;
using quantities::Square;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using testing_utilities::VanishesBefore;
using ::testing::Eq;

class SymmetricBilinearFormTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
  using Eigenworld = Frame<serialization::Frame::TestTag,
                           Inertial,
                           Handedness::Right,
                           serialization::Frame::TEST1>;

  SymmetricBilinearFormTest() {}

  template<typename Frame>
  SymmetricBilinearForm<Length, Frame, Vector> MakeSymmetricBilinearForm(
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
    return SymmetricBilinearForm<Length, Frame, Vector>(typed_matrix);
  }
};

TEST_F(SymmetricBilinearFormTest, Equality) {
  auto const f1 = MakeSymmetricBilinearForm<World>(
                      R3x3Matrix<double>({1,  2,  4},
                                         {2, -3,  5},
                                         {4,  5,  0}));
  auto const f2 = MakeSymmetricBilinearForm<World>(
                      R3x3Matrix<double>({1, 2, 3},
                                         {2, 4, 0},
                                         {3, 0, 5}));
  EXPECT_TRUE(f1 == f1);
  EXPECT_TRUE(f1 != f2);
}

TEST_F(SymmetricBilinearFormTest, UnaryOperators) {
  auto const f = MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({1,  2,  4},
                                        {2, -3,  5},
                                        {4,  5,  0}));
  EXPECT_THAT(+f,
              Eq(MakeSymmetricBilinearForm<World>(
                  R3x3Matrix<double>({1,  2,  4},
                                     {2, -3,  5},
                                     {4,  5,  0}))));
  EXPECT_THAT(-f,
              Eq(MakeSymmetricBilinearForm<World>(
                  R3x3Matrix<double>({-1, -2, -4},
                                     {-2,  3, -5},
                                     {-4, -5,  0}))));
}

TEST_F(SymmetricBilinearFormTest, BinaryOperators) {
  auto const f1 = MakeSymmetricBilinearForm<World>(
                      R3x3Matrix<double>({1,  2,  4},
                                         {2, -3,  5},
                                         {4,  5,  0}));
  auto const f2 = MakeSymmetricBilinearForm<World>(
                      R3x3Matrix<double>({1, 2, 3},
                                         {2, 4, 0},
                                         {3, 0, 5}));
  EXPECT_THAT(f1 + f2,
              Eq(MakeSymmetricBilinearForm<World>(
                  R3x3Matrix<double>({2, 4, 7},
                                     {4, 1, 5},
                                     {7, 5, 5}))));
  EXPECT_THAT(f1 - f2,
              Eq(MakeSymmetricBilinearForm<World>(
                  R3x3Matrix<double>({0,  0,  1},
                                     {0, -7,  5},
                                     {1,  5,  -5}))));
  EXPECT_THAT((-2) * f1,
              Eq(MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({-2,  -4,  -8},
                                        {-4,   6, -10},
                                        {-8, -10,   0}))));
  EXPECT_THAT(f1 * 3,
              Eq(MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({ 3,  6, 12},
                                        { 6, -9, 15},
                                        {12, 15,  0}))));
  EXPECT_THAT(f1 / 2,
              Eq(MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({0.5,    1,   2},
                                        {  1, -1.5, 2.5},
                                        {  2,  2.5,   0}))));
}

TEST_F(SymmetricBilinearFormTest, Assignment) {
  auto const f1 = MakeSymmetricBilinearForm<World>(
                      R3x3Matrix<double>({1,  2,  4},
                                         {2, -3,  5},
                                         {4,  5,  0}));
  auto const f2 = MakeSymmetricBilinearForm<World>(
                      R3x3Matrix<double>({1, 2, 3},
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
              Eq(MakeSymmetricBilinearForm<World>(
                  R3x3Matrix<double>({2, 4, 7},
                                     {4, 1, 5},
                                     {7, 5, 5}))));
  EXPECT_THAT(b,
              Eq(MakeSymmetricBilinearForm<World>(
                  R3x3Matrix<double>({0,  0,  1},
                                     {0, -7,  5},
                                     {1,  5,  -5}))));
  EXPECT_THAT(c,
              Eq(MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({ 3,  6, 12},
                                        { 6, -9, 15},
                                        {12, 15,  0}))));
  EXPECT_THAT(d,
              Eq(MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({0.5,    1,   2},
                                        {  1, -1.5, 2.5},
                                        {  2,  2.5,   0}))));
}

TEST_F(SymmetricBilinearFormTest, LinearMap) {
  auto const f = MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({1,  2,  4},
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
              Eq(MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({   2,     6,  -3.5},
                                        {   6,    18, -10.5},
                                        {-3.5, -10.5,   5}))));
}

TEST_F(SymmetricBilinearFormTest, Anticommutator) {
  auto const f = MakeSymmetricBilinearForm<World>({
      {1, 2, 4},
      {2, -3, 5},
      {4, 5, 0},
  });
  Vector<Length, World> const v({2 * Metre, 5 * Metre, 1 * Metre});
  Vector<Length, World> const w({4 * Metre, 2 * Metre, 1 * Metre});
  Bivector<Length, World> const b({1 * Metre, 3 * Metre, -5 * Metre});
  EXPECT_THAT(
      Anticommutator(f, b),
      Eq(Bivector<Square<Length>, World>(
          {11 * Pow<2>(Metre), 26 * Pow<2>(Metre), -9 * Pow<2>(Metre)})));
  EXPECT_THAT(f.Anticommutator() * b, Eq(Anticommutator(f, b)));
  EXPECT_THAT(f.Anticommutator() * b, Eq(Anticommutator(f, b)));
  EXPECT_THAT(Anticommutator(f, Wedge(v, w)),
              Eq(Wedge(f * v, w) + Wedge(v, f * w)));
}

TEST_F(SymmetricBilinearFormTest, AnticommutatorDiagonalization) {
  auto const f = MakeSymmetricBilinearForm<World>({
      {1, 0, 0},
      {0, 9, 0},
      {0, 0, 2},
  });
  using VectorEigenframe = Frame<enum class VectorTag>;
  using BivectorEigenframe = Frame<enum class BivectorTag>;
  auto const eigenvector = [](int i) {
    R3Element<double> result;
    result[i] = 1;
    return Vector<double, VectorEigenframe>(result);
  };
  auto const bieigenvector = [](int i) {
    R3Element<double> result;
    result[i] = 1;
    return Bivector<double, BivectorEigenframe>(result);
  };
  auto const vector_eigensystem = f.Diagonalize<VectorEigenframe>();
  auto const bivector_eigensystem =
      f.Anticommutator().Diagonalize<BivectorEigenframe>();
  EXPECT_THAT(vector_eigensystem.rotation(eigenvector(0)),
              Componentwise(1, 0, 0));
  EXPECT_THAT(vector_eigensystem.rotation(eigenvector(1)),
              Componentwise(0, VanishesBefore(1, 1), -1));
  EXPECT_THAT(vector_eigensystem.rotation(eigenvector(2)),
              Componentwise(0, 1, VanishesBefore(1, 1)));

  EXPECT_THAT(bivector_eigensystem.rotation(bieigenvector(0)),
              Componentwise(0, 1, 0));
  EXPECT_THAT(bivector_eigensystem.rotation(bieigenvector(1)),
              Componentwise(0, 0, -1));
  EXPECT_THAT(bivector_eigensystem.rotation(bieigenvector(2)),
              Componentwise(-1, 0, 0));
}

TEST_F(SymmetricBilinearFormTest, InnerProductForm) {
  Vector<Length, World> v1({1.0 * Metre, 3.0 * Metre, -1.0 * Metre});
  Vector<double, World> v2({2.0, 6.0, -5.0});
  auto const a = InnerProductForm<World, Vector>()(v1, v2);
  EXPECT_THAT((InnerProductForm<World, Vector>()(v1, v2)), Eq(25 * Metre));
}

TEST_F(SymmetricBilinearFormTest, Apply) {
  auto const f = MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({1,  2,  4},
                                        {2, -3,  5},
                                        {4,  5,  0}));
  Vector<Length, World> v1({1.0 * Metre, 3.0 * Metre, -1.0 * Metre});
  Vector<Length, World> v2({2.0 * Metre, 6.0 * Metre, -5.0 * Metre});
  EXPECT_THAT(f(v1, v2), Eq(-161 * Pow<3>(Metre)));
}

TEST_F(SymmetricBilinearFormTest, Serialization) {
  auto const f = MakeSymmetricBilinearForm<World>(
                     R3x3Matrix<double>({1,  2,  4},
                                        {2, -3,  5},
                                        {4,  5,  0}));
  serialization::SymmetricBilinearForm message1;
  f.WriteToMessage(&message1);
  EXPECT_TRUE(message1.has_frame());
  EXPECT_TRUE(message1.has_matrix());
  auto const g =
      SymmetricBilinearForm<Length, World, Vector>::ReadFromMessage(message1);
  EXPECT_EQ(f, g);
  serialization::SymmetricBilinearForm message2;
  g.WriteToMessage(&message2);
  EXPECT_THAT(message2, EqualsProto(message1));
}

TEST_F(SymmetricBilinearFormTest, Diagonalize) {
  // A simple test where it's clear what the eigensystem should be.
  {
    auto const f = MakeSymmetricBilinearForm<World>(
                       R3x3Matrix<double>({1,  0,  0},
                                          {0, -3,  0},
                                          {0,  0,  2}));
    auto const f_eigensystem = f.Diagonalize<Eigenworld>();

    Vector<double, Eigenworld> const e₀({1, 0, 0});
    Vector<double, Eigenworld> const e₁({0, 1, 0});
    Vector<double, Eigenworld> const e₂({0, 0, 1});
    EXPECT_THAT(f_eigensystem.form(e₀, e₀), AlmostEquals(-3 * Metre, 1));
    EXPECT_THAT(f_eigensystem.form(e₀, e₁), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₀, e₂), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₁, e₀), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₁, e₁), AlmostEquals(1 * Metre, 6));
    EXPECT_THAT(f_eigensystem.form(e₁, e₂), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₂, e₀), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₂, e₁), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₂, e₂), AlmostEquals(2 * Metre, 1));

    Vector<double, World> const w₀({ 0, 1, 0});
    Vector<double, World> const w₁({-1, 0, 0});
    Vector<double, World> const w₂({ 0, 0, 1});
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₀),
                Componentwise(AlmostEquals(1, 0),
                              VanishesBefore(1, 1),
                              VanishesBefore(1, 0)));
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₁),
                Componentwise(VanishesBefore(1, 2),
                              AlmostEquals(1, 0),
                              VanishesBefore(1, 0)));
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₂),
                Componentwise(VanishesBefore(1, 0),
                              VanishesBefore(1, 0),
                              AlmostEquals(1, 0)));
  }

  // A complex test where the eigensystem was computed using Mathematica.
  {
    auto const f = MakeSymmetricBilinearForm<World>(
                       R3x3Matrix<double>({1,  3, 2},
                                          {3, -5, 7},
                                          {2,  7, 0}));
    auto const f_eigensystem = f.Diagonalize<Eigenworld>();

    Vector<double, Eigenworld> const e₀({1, 0, 0});
    Vector<double, Eigenworld> const e₁({0, 1, 0});
    Vector<double, Eigenworld> const e₂({0, 0, 1});
    EXPECT_THAT(f_eigensystem.form(e₀, e₀),
                AlmostEquals(-10.096452436666494320 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₁, e₁),
                AlmostEquals(-0.79093267638983993780 * Metre, 34));
    EXPECT_THAT(f_eigensystem.form(e₂, e₂),
                AlmostEquals(6.8873851130563342581 * Metre, 1));

    Vector<double, World> const w₀({-0.12466193785000435776,
                                     0.82678695026555030329,
                                    -0.54852779339070148904});
    Vector<double, World> const w₁({-0.85579200995087470058,
                                     0.19015172549439114900,
                                     0.48110534916559355649});
    Vector<double, World> const w₂({0.50207513078793658603,
                                    0.52940122795673242036,
                                    0.68385298338325629274});
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₀),
                Componentwise(AlmostEquals(1, 0),
                              VanishesBefore(1, 0),
                              VanishesBefore(1, 0)));
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₁),
                Componentwise(VanishesBefore(1, 0),
                              AlmostEquals(1, 0),
                              VanishesBefore(1, 0)));
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₂),
                Componentwise(VanishesBefore(1, 1),
                              VanishesBefore(1, 2),
                              AlmostEquals(1, 0)));
  }
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
