
#include "geometry/symmetric_bilinear_form.hpp"

#include <random>

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
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace geometry {
namespace internal_symmetric_bilinear_form {

using quantities::Length;
using quantities::Pow;
using quantities::Square;
using quantities::si::Metre;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::VanishesBefore;
using testing_utilities::operator""_⑴;
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
  auto const f = MakeSymmetricBilinearForm<World>({{1,  2, 4},
                                                   {2, -3, 5},
                                                   {4,  5, 0}});
  Vector<Length, World> const v({2 * Metre, 5 * Metre, 1 * Metre});
  Vector<Length, World> const w({4 * Metre, 2 * Metre, 1 * Metre});
  Bivector<Length, World> const a({3 * Metre, 8 * Metre, 0 * Metre});
  Bivector<Length, World> const b({1 * Metre, 3 * Metre, -5 * Metre});
  EXPECT_THAT(
      Anticommutator(f, b),
      Eq(Bivector<Square<Length>, World>(
          {11 * Pow<2>(Metre), 26 * Pow<2>(Metre), -9 * Pow<2>(Metre)})));

  EXPECT_THAT(f.Anticommutator() * b, Eq(Anticommutator(f, b)));
  EXPECT_THAT(f.Anticommutator()(a, b),
              Eq(InnerProduct(a, Anticommutator(f, b))));
  EXPECT_THAT(Anticommutator(f, Wedge(v, w)),
              Eq(Wedge(f * v, w) + Wedge(v, f * w)));
  EXPECT_THAT(f.Anticommutator().AnticommutatorInverse(), Eq(f));
}

TEST_F(SymmetricBilinearFormTest, AnticommutatorDiagonalization) {
  auto const f = MakeSymmetricBilinearForm<World>({{1, 0, 0},
                                                   {0, 9, 0},
                                                   {0, 0, 2}});
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
                              VanishesBefore(1, 1)));
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₁),
                Componentwise(VanishesBefore(1, 0),
                              AlmostEquals(1, 0),
                              VanishesBefore(1, 0)));
    EXPECT_THAT(f_eigensystem.rotation.Inverse()(w₂),
                Componentwise(VanishesBefore(1, 2),
                              VanishesBefore(1, 2),
                              AlmostEquals(1, 0)));
  }

  // A degenerate case.
  {
    auto const f = MakeSymmetricBilinearForm<World>(R3x3Matrix<double>(
        {2, 0, 0},
        {0, 2, 0},
        {0, 0, 2}));
    auto const f_eigensystem = f.Diagonalize<Eigenworld>();

    Vector<double, Eigenworld> const e₀({1, 0, 0});
    Vector<double, Eigenworld> const e₁({0, 1, 0});
    Vector<double, Eigenworld> const e₂({0, 0, 1});
    EXPECT_THAT(f_eigensystem.form(e₀, e₀), AlmostEquals(2 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₀, e₁), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₀, e₂), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₁, e₀), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₁, e₁), AlmostEquals(2 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₁, e₂), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₂, e₀), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₂, e₁), AlmostEquals(0 * Metre, 0));
    EXPECT_THAT(f_eigensystem.form(e₂, e₂), AlmostEquals(2 * Metre, 0));
  }

  // Random matrices.
#if !defined(_DEBUG)
  {
    std::mt19937_64 random(42);
    std::uniform_real_distribution<> inertia_tensor_distribution(-1000.0,
                                                                 1000.0);
    for (int i = 0; i < 1'000'000; ++i) {
      double const i00 = inertia_tensor_distribution(random);
      double const i01 = inertia_tensor_distribution(random);
      double const i02 = inertia_tensor_distribution(random);
      double const i11 = inertia_tensor_distribution(random);
      double const i12 = inertia_tensor_distribution(random);
      double const i22 = inertia_tensor_distribution(random);
      auto const f = MakeSymmetricBilinearForm<World>(
          R3x3Matrix<double>({i00, i01, i02},
                             {i01, i11, i12},
                             {i02, i12, i22}));
      auto const f_eigensystem = f.Diagonalize<Eigenworld>();

      EXPECT_THAT(f_eigensystem.rotation.quaternion().Norm(),
                  AlmostEquals(1.0, 0, 3)) << f;
    }
  }
#endif

  // A poorly-conditioned case that used to yield a non-unit quaternion because
  // of non-orthogonal eigenvectors.
  {
    auto const f = MakeSymmetricBilinearForm<World>(
        R3x3Matrix<double>({{+2.86210785240963560e+02,
                             +4.96371874564241580e+02,
                             +5.70775558185911450e+02},
                            {+4.96371874564241580e+02,
                             +7.16992846801826317e+02,
                             +8.67180777241129135e+02},
                            {+5.70775558185911450e+02,
                             +8.67180777241129135e+02,
                             +9.56737022360950050e+02}}));
    auto const f_eigensystem = f.Diagonalize<Eigenworld>();

    EXPECT_THAT(f_eigensystem.rotation.quaternion().Norm(),
                AlmostEquals(1.0, 1)) << f;
  }

  // A case with two eigenvalues that are very close (exact relative error
  // 2e-15) that used to yield a non-unit quaternion because of a missing
  // normalization.  Found in game in #2611.
  {
    auto const f = MakeSymmetricBilinearForm<World>(
        R3x3Matrix<double>({{+6.25360854308065672e+00,
                             +2.24243333089292812e-01,
                             +1.68316543009972008e-02},
                            {+2.24243333089292812e-01,
                             +1.09414207983843497e+01,
                             +3.52669451554594282e-01},
                            {+1.68316543009972008e-02,
                             +3.52669451554594282e-01,
                             +6.26937749903824937e+00}}));
    auto const f_eigensystem = f.Diagonalize<Eigenworld>();

    EXPECT_THAT(f_eigensystem.rotation.quaternion().Norm(),
                AlmostEquals(1.0, 0)) << f;
    Vector<double, Eigenworld> const e₀({1, 0, 0});
    Vector<double, Eigenworld> const e₁({0, 1, 0});
    Vector<double, Eigenworld> const e₂({0, 0, 1});

    // Real eigenvectors obtained with Mathematica.  Note how bad the first two
    // eigenvectors are: that's because the object is very nearly a disc, and
    // our computation yields eigenvectors that are roughly 45° from the real
    // ones.  But at least they are in the right plane and determine a rotation
    // that correctly aligns the isolated eigenvector.
    EXPECT_THAT(f_eigensystem.rotation(e₀),
                Componentwise(
                    AbsoluteErrorFrom(-0.71267592684216601514, IsNear(0.7_⑴)),
                    AbsoluteErrorFrom(+0.086267747252993862323, IsNear(0.01_⑴)),
                    AbsoluteErrorFrom(-0.69616872888945048034, IsNear(0.3_⑴))));
    EXPECT_THAT(f_eigensystem.rotation(e₁),
                Componentwise(
                    AbsoluteErrorFrom(-0.69988076924532898781, IsNear(0.3_⑴)),
                    AbsoluteErrorFrom(-0.02018794239465783510, IsNear(0.07_⑴)),
                    AbsoluteErrorFrom(+0.71397433835008141314, IsNear(0.7_⑴))));
    EXPECT_THAT(f_eigensystem.rotation(e₂),
                Componentwise(AlmostEquals(+0.04753874357012595212, 10),
                              AlmostEquals(+0.99606742882485799451, 1),
                              AlmostEquals(+0.07476459786563599011, 4)));
  }

  // A case similar to the previous one, but constructed so that the two largest
  // eigenvalues are very close (a needle).
  {
    auto const f = MakeSymmetricBilinearForm<World>(
        R3x3Matrix<double>({{+3.02958892130780040082,
                             +0.34454179629510833170,
                             -3.74544524094010311734},
                            {+0.34454179629510833170,
                             +9.98296957696546981827,
                             +0.18513433665111470179},
                            {-3.74544524094010311734,
                             +0.18513433665111470179,
                             +7.98744150172682978089}}));
    auto const f_eigensystem = f.Diagonalize<Eigenworld>();

    EXPECT_THAT(f_eigensystem.rotation.quaternion().Norm(),
                AlmostEquals(1.0, 0)) << f;
    Vector<double, Eigenworld> const e₀({1, 0, 0});
    Vector<double, Eigenworld> const e₁({0, 1, 0});
    Vector<double, Eigenworld> const e₂({0, 0, 1});

    // Same comment as above regarding the last two eigenvectors.
    EXPECT_THAT(f_eigensystem.rotation(e₀),
                Componentwise(AlmostEquals(+0.88005120297326585654, 1),
                              AlmostEquals(-0.04350022098854655144, 1),
                              AlmostEquals(+0.47288223789782508101, 2)));
    EXPECT_THAT(f_eigensystem.rotation(e₁),
                Componentwise(
                    AbsoluteErrorFrom(-0.42509236497268917818, IsNear(0.07_⑴)),
                    AbsoluteErrorFrom(+0.37170808088366222608, IsNear(1.1_⑴)),
                    AbsoluteErrorFrom(+0.82530575173550731715, IsNear(0.2_⑴))));
    EXPECT_THAT(f_eigensystem.rotation(e₂),
                Componentwise(
                    AbsoluteErrorFrom(+0.21167513171658507292, IsNear(0.1_⑴)),
                    AbsoluteErrorFrom(+0.92732994849715299942, IsNear(1.6_⑴)),
                    AbsoluteErrorFrom(-0.30863053191969508327, IsNear(0.3_⑴))));
  }
}

}  // namespace internal_symmetric_bilinear_form
}  // namespace geometry
}  // namespace principia
