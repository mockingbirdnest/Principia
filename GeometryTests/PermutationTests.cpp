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
  struct World;
  typedef Permutation<quantities::Length, World, World> P;
 public:
  TEST_METHOD_INITIALIZE(Initialize) {
    vector_.reset(new Vector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
  }

  TEST_METHOD(XYZ) {
    Vector<quantities::Length, World> const permuted = 
        P(P::XYZ) * vector_->coordinates;
    AssertEqual(1.0 * Metre, permuted.coordinates.x);
    AssertEqual(2.0 * Metre, permuted.coordinates.y);
    AssertEqual(3.0 * Metre, permuted.coordinates.z);
  }

  TEST_METHOD(YZX) {
    Vector<quantities::Length, World> const permuted = 
        P(P::YZX) * vector_->coordinates;
    AssertEqual(2.0 * Metre, permuted.coordinates.x);
    AssertEqual(3.0 * Metre, permuted.coordinates.y);
    AssertEqual(1.0 * Metre, permuted.coordinates.z);
  }

  TEST_METHOD(ZXY) {
    Vector<quantities::Length, World> const permuted = 
        P(P::ZXY) * vector_->coordinates;
    AssertEqual(3.0 * Metre, permuted.coordinates.x);
    AssertEqual(1.0 * Metre, permuted.coordinates.y);
    AssertEqual(2.0 * Metre, permuted.coordinates.z);
  }

  TEST_METHOD(XZY) {
    Vector<quantities::Length, World> const permuted = 
        P(P::XZY) * vector_->coordinates;
    AssertEqual(1.0 * Metre, permuted.coordinates.x);
    AssertEqual(3.0 * Metre, permuted.coordinates.y);
    AssertEqual(2.0 * Metre, permuted.coordinates.z);
  }

  TEST_METHOD(ZYX) {
    Vector<quantities::Length, World> const permuted = 
        P(P::ZYX) * vector_->coordinates;
    AssertEqual(3.0 * Metre, permuted.coordinates.x);
    AssertEqual(2.0 * Metre, permuted.coordinates.y);
    AssertEqual(1.0 * Metre, permuted.coordinates.z);
  }

  TEST_METHOD(YXZ) {
    Vector<quantities::Length, World> const permuted = 
        P(P::YXZ) * vector_->coordinates;
    AssertEqual(2.0 * Metre, permuted.coordinates.x);
    AssertEqual(1.0 * Metre, permuted.coordinates.y);
    AssertEqual(3.0 * Metre, permuted.coordinates.z);
  }

 private:
  std::unique_ptr<Vector<quantities::Length, World>> vector_;
};

}  // namespace geometry
}  // namespace principia
