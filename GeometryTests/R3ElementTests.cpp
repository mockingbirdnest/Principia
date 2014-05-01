
#include "stdafx.hpp"

#include <functional>

#include "CppUnitTest.h"
#include "..\TestUtilities\TestUtilities.hpp"
#include "..\TestUtilities\QuantityComparisons.hpp"
#include "..\TestUtilities\GeometryComparisons.hpp"
#include "..\TestUtilities\Algebra.hpp"

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"
#include "..\Quantities\SI.hpp"
#include "..\Quantities\Constants.hpp"
#include "..\Quantities\UK.hpp"
#include "..\Quantities\BIPM.hpp"
#include "..\Quantities\Astronomy.hpp"
#include "..\Geometry\R3Element.hpp"
#include "..\Geometry\Grassmann.hpp"
#include "..\Quantities\ElementaryFunctions.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Principia {
namespace Geometry {
namespace {

using namespace TestUtilities;

using namespace Quantities;
using namespace SI;
using namespace UK;
using namespace BIPM;
using namespace Astronomy;
using namespace Constants;

TEST_CLASS(R3ElementTests) {
 public:
  TEST_METHOD(Dumb3Vector) {
    R3Element<Speed> nullVector(0 * Metre / Second,
                                0 * Metre / Second,
                                0 * Metre / Second);
    R3Element<Speed> u(1 * Metre / Second,
                       120 * Kilo(Metre) / Hour,
                       -SpeedOfLight);
    R3Element<Speed> v(-20 * Knot,
                       2 * π * AstronomicalUnit / JulianYear,
                       1 * Admiralty::NauticalMile / Hour);
    R3Element<Speed> w(-1 * Mile / Hour, -2 * Foot / Second, -3 * Knot);
    R3Element<Speed> a(88 * Mile / Hour, 300 * Metre / Second, 46 * Knot);
    AssertEqual((e * Dimensionless(42)) * v, e * (Dimensionless(42) * v));
    TestVectorSpace<R3Element<Speed>,
                    Dimensionless>(nullVector, u, v, w, Dimensionless(0),
                                   Dimensionless(1), e, Dimensionless(42));
    TestAlternatingBilinearMap(Cross<Speed, Speed>, u,
                               v, w, a, Dimensionless(42));
    TestSymmetricPositiveDefiniteBilinearMap(Dot<Speed, Speed>,
                                             u, v, w, a, Dimensionless(42));
  }
};

}  // namepsace
}  // namespace Geometry
}  // namespace Principia
