#include "stdafx.hpp"

#include <CppUnitTest.h>
#include <memory>

#include "Geometry/Grassmann.hpp"
#include "Geometry/Permutation.hpp"
#include "Quantities/Quantities.hpp"
#include "Quantities/SI.hpp"
#include "TestUtilities/TestUtilities.hpp"

namespace principia {
namespace geometry {

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace si;
using namespace test_utilities;

TEST_CLASS(PermutationTests) {
 private:
  struct World;
  typedef Permutation<quantities::Length, World, World> P;

  std::unique_ptr<Vector<quantities::Length, World>> vector_;
  std::unique_ptr<Bivector<quantities::Length, World>> bivector_;
  std::unique_ptr<Trivector<quantities::Length, World>> trivector_;

 public:
  TEST_METHOD_INITIALIZE(Initialize) {
    vector_.reset(new Vector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    bivector_.reset(new Bivector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    trivector_.reset(new Trivector<quantities::Length, World>(4.0 * Metre));
  }

  TEST_METHOD(XYZ) {
    R3Element<quantities::Length> const permuted = 
        P(P::XYZ) * vector_->coordinates;
    AssertEqual(1.0 * Metre, permuted.x);
    AssertEqual(2.0 * Metre, permuted.y);
    AssertEqual(3.0 * Metre, permuted.z);
  }

  TEST_METHOD(YZX) {
    R3Element<quantities::Length> const permuted = 
        P(P::YZX) * vector_->coordinates;
    AssertEqual(2.0 * Metre, permuted.x);
    AssertEqual(3.0 * Metre, permuted.y);
    AssertEqual(1.0 * Metre, permuted.z);
  }

  TEST_METHOD(ZXY) {
    R3Element<quantities::Length> const permuted = 
        P(P::ZXY) * vector_->coordinates;
    AssertEqual(3.0 * Metre, permuted.x);
    AssertEqual(1.0 * Metre, permuted.y);
    AssertEqual(2.0 * Metre, permuted.z);
  }

  TEST_METHOD(XZY) {
    R3Element<quantities::Length> const permuted = 
        P(P::XZY) * vector_->coordinates;
    AssertEqual(1.0 * Metre, permuted.x);
    AssertEqual(3.0 * Metre, permuted.y);
    AssertEqual(2.0 * Metre, permuted.z);
  }

  TEST_METHOD(ZYX) {
    R3Element<quantities::Length> const permuted = 
        P(P::ZYX) * vector_->coordinates;
    AssertEqual(3.0 * Metre, permuted.x);
    AssertEqual(2.0 * Metre, permuted.y);
    AssertEqual(1.0 * Metre, permuted.z);
  }

  TEST_METHOD(YXZ) {
    R3Element<quantities::Length> const permuted = 
        P(P::YXZ) * vector_->coordinates;
    AssertEqual(2.0 * Metre, permuted.x);
    AssertEqual(1.0 * Metre, permuted.y);
    AssertEqual(3.0 * Metre, permuted.z);
  }

  TEST_METHOD(Determinant) {
    P xyz(P::XYZ);
    P yzx(P::YZX);
    P zxy(P::ZXY);
    P xzy(P::XZY);
    P zyx(P::ZYX);
    P yxz(P::YXZ);
    AssertTrue(xyz.Determinant().Positive());
    AssertTrue(yzx.Determinant().Positive());
    AssertTrue(zxy.Determinant().Positive());
    AssertTrue(xzy.Determinant().Negative());
    AssertTrue(zyx.Determinant().Negative());
    AssertTrue(yxz.Determinant().Negative());
  }

  TEST_METHOD(AppliedToVector) {
    Vector<quantities::Length, World> permuted_yzx =
        P(P::YZX)(*vector_);
    AssertEqual(2.0 * Metre, permuted_yzx.coordinates.x);
    AssertEqual(3.0 * Metre, permuted_yzx.coordinates.y);
    AssertEqual(1.0 * Metre, permuted_yzx.coordinates.z);
    Vector<quantities::Length, World> permuted_yxz =
        P(P::YXZ)(*vector_);
    AssertEqual(2.0 * Metre, permuted_yxz.coordinates.x);
    AssertEqual(1.0 * Metre, permuted_yxz.coordinates.y);
    AssertEqual(3.0 * Metre, permuted_yxz.coordinates.z);
  }

  TEST_METHOD(AppliedToBivector) {
    Bivector<quantities::Length, World> permuted_zxy =
        P(P::ZXY)(*bivector_);
    AssertEqual(3.0 * Metre, permuted_zxy.coordinates.x);
    AssertEqual(1.0 * Metre, permuted_zxy.coordinates.y);
    AssertEqual(2.0 * Metre, permuted_zxy.coordinates.z);
    Bivector<quantities::Length, World> permuted_zyx =
        P(P::ZYX)(*bivector_);
    AssertEqual(-3.0 * Metre, permuted_zyx.coordinates.x);
    AssertEqual(-2.0 * Metre, permuted_zyx.coordinates.y);
    AssertEqual(-1.0 * Metre, permuted_zyx.coordinates.z);
  }

  TEST_METHOD(AppliedToTrivector) {
    Trivector<quantities::Length, World> permuted_xyz =
        P(P::XYZ)(*trivector_);
    AssertEqual(4.0 * Metre, permuted_xyz.coordinates);
    Trivector<quantities::Length, World> permuted_xzy =
        P(P::XZY)(*trivector_);
    AssertEqual(-4.0 * Metre, permuted_xzy.coordinates);
  }
};

}  // namespace geometry
}  // namespace principia
