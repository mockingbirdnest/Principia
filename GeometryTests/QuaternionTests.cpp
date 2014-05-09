#include "stdafx.hpp"

#include <CppUnitTest.h>

#include "Geometry/Quaternion.hpp"
#include "TestUtilities/TestUtilities.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace principia {
namespace geometry {

using namespace test_utilities;

TEST_CLASS(QuaternionTests) {
 public:
  struct World;
  typedef Permutation<World, World> P;

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

  TEST_METHOD(Integer) {
    Sign const positive(1);
    Sign const negative(-1);
    AssertTrue(positive.Positive());
    AssertFalse(positive.Negative());
    AssertFalse(negative.Positive());
    AssertTrue(negative.Negative());
  }

  TEST_METHOD(SignMultiplication) {
    Sign const positive(1);
    Sign const negative(-1);
    AssertTrue((positive * positive).Positive());
    AssertTrue((positive * negative).Negative());
    AssertTrue((negative * positive).Negative());
    AssertTrue((negative * negative).Positive());
  }

  TEST_METHOD(ScalarMultiplication) {
    Sign const positive(1);
    Sign const negative(-1);
    AssertEqual(3, positive * 3);
    AssertEqual(-3, positive * -3);
    AssertEqual(-3, negative * 3);
    AssertEqual(3, negative * -3);
  }
};

}  // namespace geometry
}  // namespace principia
