#include "stdafx.hpp"
#include <float.h>

#include <CppUnitTest.h>
#include <memory>

#include "Geometry/Grassmann.hpp"
#include "Geometry/Quaternion.hpp"
#include "Geometry/Rotation.hpp"
#include "Quantities/Quantities.hpp"
#include "Quantities/SI.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace geometry {

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace si;
using namespace test_utilities;

TEST_CLASS(RotationTests) {
 private:
  struct World;
  typedef Rotation<World, World> R;

  std::unique_ptr<Vector<quantities::Length, World>> vector_;
  std::unique_ptr<Bivector<quantities::Length, World>> bivector_;
  std::unique_ptr<Trivector<quantities::Length, World>> trivector_;
  std::unique_ptr<R> rotation_a_;
  std::unique_ptr<R> rotation_b_;
  std::unique_ptr<R> rotation_c_;

 public:
  TEST_METHOD_INITIALIZE(Initialize) {
    vector_.reset(new Vector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    bivector_.reset(new Bivector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    trivector_.reset(new Trivector<quantities::Length, World>(4.0 * Metre));
    rotation_a_.reset(
      new R(120 * si::Degree, 
            Vector<quantities::Dimensionless, World>({1, 1, 1})));
    rotation_b_.reset(
      new R(90 * si::Degree, 
            Vector<quantities::Dimensionless, World>({1, 0, 0})));
    rotation_c_.reset(
      new R(30 * si::Degree, 
            Vector<quantities::Dimensionless, World>({1, 1, 1})));
  }

  TEST_METHOD(Identity) {
    R3Element<quantities::Length> const rotated =
        R::Identity()(vector_->coordinates());
    AssertEqual(1.0 * Metre, rotated.x);
    AssertEqual(2.0 * Metre, rotated.y);
    AssertEqual(3.0 * Metre, rotated.z);
  }

  TEST_METHOD(AppliedToVector) {
    Vector<quantities::Length, World> const rotated = (*rotation_a_)(*vector_);
    AssertEqual(3.0 * Metre, rotated.coordinates().x);
    AssertEqualWithin(1.0 * Metre, rotated.coordinates().y, 3 * DBL_EPSILON);
    AssertEqual(2.0 * Metre, rotated.coordinates().z);
  }

  TEST_METHOD(AppliedToBivector) {
    Bivector<quantities::Length, World> const rotated =
        (*rotation_a_)(*bivector_);
    AssertEqual(3.0 * Metre, rotated.coordinates().x);
    AssertEqualWithin(1.0 * Metre, rotated.coordinates().y, 3 * DBL_EPSILON);
    AssertEqual(2.0 * Metre, rotated.coordinates().z);
  }

  TEST_METHOD(AppliedToTrivector) {
    Trivector<quantities::Length, World> const rotated =
        (*rotation_a_)(*trivector_);
    AssertEqual(4.0 * Metre, rotated.coordinates());
  }

  TEST_METHOD(Determinant) {
    Sign const determinant = rotation_a_->Determinant();
    AssertTrue(determinant.Positive());
  }

  TEST_METHOD(Inverse) {
    R const inverse = rotation_a_->Inverse();
    Vector<quantities::Length, World> const rotated = inverse(*vector_);
    AssertEqualWithin(2.0 * Metre, rotated.coordinates().x, 3 * DBL_EPSILON);
    AssertEqualWithin(3.0 * Metre, rotated.coordinates().y, 3 * DBL_EPSILON);
    AssertEqual(1.0 * Metre, rotated.coordinates().z);
  }

  TEST_METHOD(Composition) {
    R const rotation_ab_ = *rotation_a_ * *rotation_b_;
    Vector<quantities::Length, World> const rotated = rotation_ab_(*vector_);
    AssertEqualWithin(2.0 * Metre, rotated.coordinates().x, 3 * DBL_EPSILON);
    AssertEqualWithin(1.0 * Metre, rotated.coordinates().y, 3 * DBL_EPSILON);
    AssertEqual(-3.0 * Metre, rotated.coordinates().z);
  }
};

}  // namespace geometry
}  // namespace principia
