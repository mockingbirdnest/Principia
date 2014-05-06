#include "stdafx.hpp"

#include <CppUnitTest.h>

#include "Geometry/Grassmann.hpp"
#include "Geometry/R3Element.hpp"
#include "Quantities/Astronomy.hpp"
#include "Quantities/BIPM.hpp"
#include "Quantities/Constants.hpp"
#include "Quantities/Dimensionless.hpp"
#include "Quantities/ElementaryFunctions.hpp"
#include "Quantities/Quantities.hpp"
#include "Quantities/SI.hpp"
#include "Quantities/UK.hpp"
#include "TestUtilities/Algebra.hpp"
#include "TestUtilities/GeometryComparisons.hpp"
#include "TestUtilities/QuantityComparisons.hpp"
#include "TestUtilities/TestUtilities.hpp"


namespace principia {
namespace geometry {

using namespace astronomy;
using namespace bipm;
using namespace constants;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace quantities;
using namespace si;
using namespace test_utilities;
using namespace uk;

TEST_CLASS(R3ElementTests) {
  R3Element<Speed> const null_velocity_ =
      R3Element<Speed>(0 * Knot, 0 * Knot, 0 * Knot);
  R3Element<Speed> const u_ =
      R3Element<Speed>(3 * Knot, -42 * Parsec / JulianYear, 0 * Knot);
  R3Element<Speed> const v_ =
      R3Element<Speed>(-π * SpeedOfLight, -e * Kilo(Metre) / Hour, -1 * Knot);
  R3Element<Speed> const w_ =
      R3Element<Speed>(2 * Mile / Hour, 2 * Furlong / Day, 2 * Rod / Minute);
  R3Element<Speed> const a_ =
      R3Element<Speed>(88 * Mile / Hour, 300 * Metre / Second, 46 * Knot);

 public:
  TEST_METHOD(Dumb3Vector) {
    AssertEqual((e * Dimensionless(42)) * v_, e * (Dimensionless(42) * v_));
    TestVectorSpace<R3Element<Speed>, Dimensionless>(null_velocity_, u_, v_,
                                                     w_, Dimensionless(0),
                                                     Dimensionless(1), e,
                                                     Dimensionless(42));
    TestAlternatingBilinearMap(Cross<Speed, Speed>, u_, v_, w_, a_,
                               Dimensionless(42));
    TestSymmetricPositiveDefiniteBilinearMap(Dot<Speed, Speed>, u_, v_, w_, a_,
                                             Dimensionless(42));
  }

  TEST_METHOD(MixedProduct) {
    auto left_time_multiplication = [](Time left, R3Element<Speed> right) {
      return left * right;
    };
    auto right_time_multiplication = [](R3Element<Speed> left, Time right) {
      return left * right;
    };
    TestBilinearMap(left_time_multiplication, 1 * Second, 1 * JulianYear, u_,
                    v_, Dimensionless(42));
    TestBilinearMap(right_time_multiplication, w_, a_, -1 * Day,
                    1 * Parsec / SpeedOfLight, Dimensionless(-π));
    Time t = -3 * Second;
    AssertEqual(t * u_, u_ * t);
    AssertEqual((u_ * t) / t, u_);
  }
};

}  // namespace geometry
}  // namespace principia
