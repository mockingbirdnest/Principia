#include "stdafx.hpp"

#include <CppUnitTest.h>
#include <memory>

#include "geometry/quaternion.hpp"
#include "TestUtilities/GeometryComparisons.hpp"
#include "TestUtilities/TestUtilities.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace principia {
namespace geometry {

using namespace test_utilities;

TEST_CLASS(QuaternionTests) {
 private:
  std::unique_ptr<Quaternion> q1_;
  std::unique_ptr<Quaternion> q2_;
  std::unique_ptr<Quaternion> q3_;

 public:
  TEST_METHOD_INITIALIZE(Initialize) {
    q1_.reset(new Quaternion(1, {1, -1, -1}));
    q2_.reset(new Quaternion(-2, {1, -3, 4}));
    q3_.reset(new Quaternion(-8));
  }

  TEST_METHOD(RealPart) {
    AssertEqual<quantities::Dimensionless>(1, q1_->real_part());
    AssertEqual<quantities::Dimensionless>(-2, q2_->real_part());
    AssertEqual<quantities::Dimensionless>(-8, q3_->real_part());
  }

  TEST_METHOD(ImaginaryPart) {
    AssertEqual<quantities::Dimensionless>({1, -1, -1}, q1_->imaginary_part());
    AssertEqual<quantities::Dimensionless>({1, -3, 4}, q2_->imaginary_part());
    AssertEqual<quantities::Dimensionless>({0, 0, 0}, q3_->imaginary_part());
  }

  TEST_METHOD(Conjugate) {
    Quaternion const p1(q1_->Conjugate());
    Quaternion const p2(q2_->Conjugate());
    Quaternion const p3(q3_->Conjugate());
    AssertEqual<quantities::Dimensionless>(1, p1.real_part());
    AssertEqual<quantities::Dimensionless>({-1, 1, 1}, p1.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-2, p2.real_part());
    AssertEqual<quantities::Dimensionless>({-1, 3, -4}, p2.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-8, p3.real_part());
    AssertEqual<quantities::Dimensionless>({0, 0, 0}, p3.imaginary_part());
  }

  TEST_METHOD(Inverse) {
    Quaternion const p1(q1_->Inverse());
    Quaternion const p2(q2_->Inverse());
    Quaternion const p3(q3_->Inverse());
    Quaternion const r1(p1 * *q1_);
    AssertEqual<quantities::Dimensionless>(0.25, p1.real_part());
    AssertEqual<quantities::Dimensionless>({-0.25, 0.25, 0.25},
                                           p1.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-1.0 / 15.0, p2.real_part());
    AssertEqual<quantities::Dimensionless>({-1.0 / 30.0, 0.1, -2.0 / 15.0},
                                           p2.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-0.125, p3.real_part());
    AssertEqual<quantities::Dimensionless>({0, 0, 0}, p3.imaginary_part());

    AssertEqual<quantities::Dimensionless>(1, r1.real_part());
    AssertEqual<quantities::Dimensionless>({0, 0, 0}, r1.imaginary_part());
  }

  TEST_METHOD(UnaryOperators) {
    Quaternion const a(+*q2_);
    Quaternion const b(-*q2_);
    AssertEqual<quantities::Dimensionless>(-2, a.real_part());
    AssertEqual<quantities::Dimensionless>({1, -3, 4}, a.imaginary_part());
    AssertEqual<quantities::Dimensionless>(2, b.real_part());
    AssertEqual<quantities::Dimensionless>({-1, 3, -4}, b.imaginary_part());
  }

  TEST_METHOD(BinaryOperators) {
    Quaternion const a(*q2_ + *q1_);
    Quaternion const b(*q2_ - *q1_);
    Quaternion const c(*q2_ * *q1_);
    Quaternion const d(*q2_ / *q1_);
    AssertEqual<quantities::Dimensionless>(-1, a.real_part());
    AssertEqual<quantities::Dimensionless>({2, -4, 3}, a.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-3, b.real_part());
    AssertEqual<quantities::Dimensionless>({0, -2, 5}, b.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-2, c.real_part());
    AssertEqual<quantities::Dimensionless>({6, 4, 8}, c.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-0.5, d.real_part());
    AssertEqual<quantities::Dimensionless>({-1, -2.5, 0}, d.imaginary_part());
  }

  TEST_METHOD(ScalarMultiplicationDivision) {
    Quaternion const a(3 * *q2_);
    Quaternion const b(*q2_ * -5);
    Quaternion const c(*q2_ / 4);
    AssertEqual<quantities::Dimensionless>(-6, a.real_part());
    AssertEqual<quantities::Dimensionless>({3, -9, 12}, a.imaginary_part());
    AssertEqual<quantities::Dimensionless>(10, b.real_part());
    AssertEqual<quantities::Dimensionless>({-5, 15, -20}, b.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-0.5, c.real_part());
    AssertEqual<quantities::Dimensionless>({0.25, -0.75, 1}, c.imaginary_part());
  }

  TEST_METHOD(Assignment) {
    Quaternion a = *q2_;
    Quaternion b = *q2_;
    Quaternion c = *q2_;
    Quaternion d = *q2_;
    Quaternion e = *q2_;
    Quaternion f = *q2_;
    a += *q1_;
    b -= *q1_;
    c *= *q1_;
    d /= *q1_;
    e *= 3;
    f /= 4;
    AssertEqual<quantities::Dimensionless>(-1, a.real_part());
    AssertEqual<quantities::Dimensionless>({2, -4, 3}, a.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-3, b.real_part());
    AssertEqual<quantities::Dimensionless>({0, -2, 5}, b.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-2, c.real_part());
    AssertEqual<quantities::Dimensionless>({6, 4, 8}, c.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-0.5, d.real_part());
    AssertEqual<quantities::Dimensionless>({-1, -2.5, 0}, d.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-6, e.real_part());
    AssertEqual<quantities::Dimensionless>({3, -9, 12}, e.imaginary_part());
    AssertEqual<quantities::Dimensionless>(-0.5, f.real_part());
    AssertEqual<quantities::Dimensionless>({0.25, -0.75, 1}, f.imaginary_part());
  }

};

}  // namespace geometry
}  // namespace principia
