
#include "geometry/quaternion.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/algebra.hpp"

namespace principia {

using testing_utilities::TestSkewField;
using ::testing::Eq;

namespace geometry {

class QuaternionTest : public testing::Test {
 protected:
  using R3 = R3Element<double>;

  QuaternionTest() {
    q1_ = Quaternion(1, {1, -1, -1});
    q2_ = Quaternion(-2, {1, -3, 4});
    q3_ = Quaternion(-8);
  }

  Quaternion q1_;
  Quaternion q2_;
  Quaternion q3_;
};

TEST_F(QuaternionTest, RealPart) {
  EXPECT_THAT(q1_.real_part(), Eq(1));
  EXPECT_THAT(q2_.real_part(), Eq(-2));
  EXPECT_THAT(q3_.real_part(), Eq(-8));
}

TEST_F(QuaternionTest, ImaginaryPart) {
  EXPECT_THAT(q1_.imaginary_part(), Eq<R3>({1, -1, -1}));
  EXPECT_THAT(q2_.imaginary_part(), Eq<R3>({1, -3, 4}));
  EXPECT_THAT(q3_.imaginary_part(), Eq<R3>({0, 0, 0}));
}

TEST_F(QuaternionTest, Conjugate) {
  Quaternion const p1 = q1_.Conjugate();
  Quaternion const p2 = q2_.Conjugate();
  Quaternion const p3 = q3_.Conjugate();
  EXPECT_THAT(p1.real_part(), Eq(1));
  EXPECT_THAT(p1.imaginary_part(), Eq<R3>({-1, 1, 1}));
  EXPECT_THAT(p2.real_part(), Eq(-2));
  EXPECT_THAT(p2.imaginary_part(), Eq<R3>({-1, 3, -4}));
  EXPECT_THAT(p3.real_part(), Eq(-8));
  EXPECT_THAT(p3.imaginary_part(), Eq<R3>({0, 0, 0}));
}

TEST_F(QuaternionTest, Inverse) {
  Quaternion const p1 = q1_.Inverse();
  Quaternion const p2 = q2_.Inverse();
  Quaternion const p3 = q3_.Inverse();
  Quaternion const r1 = p1 * q1_;
  EXPECT_THAT(p1.real_part(), Eq(0.25));
  EXPECT_THAT(p1.imaginary_part(), Eq<R3>({-0.25, 0.25, 0.25}));
  EXPECT_THAT(p2.real_part(), Eq(-1.0 / 15.0));
  EXPECT_THAT(p2.imaginary_part(), Eq<R3>({-1.0 / 30.0, 0.1, -2.0 / 15.0}));
  EXPECT_THAT(p3.real_part(), Eq(-0.125));
  EXPECT_THAT(p3.imaginary_part(), Eq<R3>({0, 0, 0}));

  EXPECT_THAT(r1.real_part(), Eq(1));
  EXPECT_THAT(r1.imaginary_part(), Eq<R3>({0, 0, 0}));
}

TEST_F(QuaternionTest, UnaryOperators) {
  Quaternion const a = +q2_;
  Quaternion const b = -q2_;
  EXPECT_THAT(a.real_part(), Eq(-2));
  EXPECT_THAT(a.imaginary_part(), Eq<R3>({1, -3, 4}));
  EXPECT_THAT(b.real_part(), Eq(2));
  EXPECT_THAT(b.imaginary_part(), Eq<R3>({-1, 3, -4}));
}

TEST_F(QuaternionTest, BinaryOperators) {
  Quaternion const a = q2_ + q1_;
  Quaternion const b = q2_ - q1_;
  Quaternion const c = q2_ * q1_;
  Quaternion const d = q2_ / q1_;
  EXPECT_THAT(a.real_part(), Eq(-1));
  EXPECT_THAT(a.imaginary_part(), Eq<R3>({2, -4, 3}));
  EXPECT_THAT(b.real_part(), Eq(-3));
  EXPECT_THAT(b.imaginary_part(), Eq<R3>({0, -2, 5}));
  EXPECT_THAT(c.real_part(), Eq(-2));
  EXPECT_THAT(c.imaginary_part(), Eq<R3>({6, 4, 8}));
  EXPECT_THAT(d.real_part(), Eq(-0.5));
  EXPECT_THAT(d.imaginary_part(), Eq<R3>({-1, -2.5, 0}));
}

TEST_F(QuaternionTest, ScalarMultiplicationDivision) {
  Quaternion const a = 3 * q2_;
  Quaternion const b = q2_ * -5;
  Quaternion const c = q2_ / 4;
  EXPECT_THAT(a.real_part(), Eq(-6));
  EXPECT_THAT(a.imaginary_part(), Eq<R3>({3, -9, 12}));
  EXPECT_THAT(b.real_part(), Eq(10));
  EXPECT_THAT(b.imaginary_part(), Eq<R3>({-5, 15, -20}));
  EXPECT_THAT(c.real_part(), Eq(-0.5));
  EXPECT_THAT(c.imaginary_part(), Eq<R3>({0.25, -0.75, 1}));
}

TEST_F(QuaternionTest, Assignment) {
  Quaternion a = q2_;
  Quaternion b = q2_;
  Quaternion c = q2_;
  Quaternion d = q2_;
  Quaternion e = q2_;
  Quaternion f = q2_;
  a += q1_;
  b -= q1_;
  c *= q1_;
  d /= q1_;
  e *= 3;
  f /= 4;
  EXPECT_THAT(a.real_part(), Eq(-1));
  EXPECT_THAT(a.imaginary_part(), Eq<R3>({2, -4, 3}));
  EXPECT_THAT(b.real_part(), Eq(-3));
  EXPECT_THAT(b.imaginary_part(), Eq<R3>({0, -2, 5}));
  EXPECT_THAT(c.real_part(), Eq(-2));
  EXPECT_THAT(c.imaginary_part(), Eq<R3>({6, 4, 8}));
  EXPECT_THAT(d.real_part(), Eq(-0.5));
  EXPECT_THAT(d.imaginary_part(), Eq<R3>({-1, -2.5, 0}));
  EXPECT_THAT(e.real_part(), Eq(-6));
  EXPECT_THAT(e.imaginary_part(), Eq<R3>({3, -9, 12}));
  EXPECT_THAT(f.real_part(), Eq(-0.5));
  EXPECT_THAT(f.imaginary_part(), Eq<R3>({0.25, -0.75, 1}));
}

TEST_F(QuaternionTest, SkewField) {
  TestSkewField<Quaternion>(
      Quaternion(0), Quaternion(1),
      q1_, q2_, q3_,
      Quaternion(6, {1, -3, 4}), Quaternion(0, {8, 9, -1}), 0, 1);
}

TEST_F(QuaternionTest, Serialization) {
  serialization::Quaternion message;
  q2_.WriteToMessage(&message);
  EXPECT_EQ(-2.0, message.real_part());
  EXPECT_TRUE(message.imaginary_part().x().has_double_());
  EXPECT_EQ(1.0, message.imaginary_part().x().double_());
  EXPECT_TRUE(message.imaginary_part().y().has_double_());
  EXPECT_EQ(-3.0, message.imaginary_part().y().double_());
  EXPECT_TRUE(message.imaginary_part().z().has_double_());
  EXPECT_EQ(4.0, message.imaginary_part().z().double_());
  Quaternion const q = Quaternion::ReadFromMessage(message);
  EXPECT_EQ(q2_, q);
}

TEST_F(QuaternionTest, Output) {
  std::cout << q2_ << "\n";
}

}  // namespace geometry
}  // namespace principia
