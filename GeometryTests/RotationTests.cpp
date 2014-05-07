#include "stdafx.hpp"

#include <CppUnitTest.h>
#include <memory>

#include "Geometry/Grassmann.hpp"
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

 public:
  TEST_METHOD_INITIALIZE(Initialize) {
    vector_.reset(new Vector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    bivector_.reset(new Bivector<quantities::Length, World>(
      R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre)));
    trivector_.reset(new Trivector<quantities::Length, World>(4.0 * Metre));
  }

  TEST_METHOD(Identity) {
    R3Element<quantities::Length> const rotated =
        R::Identity()(vector_->coordinates());
    AssertEqual(1.0 * Metre, rotated.x);
    AssertEqual(2.0 * Metre, rotated.y);
    AssertEqual(3.0 * Metre, rotated.z);
  }

};

}  // namespace geometry
}  // namespace principia
