#include "stdafx.hpp"

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
namespace GeometryTests
{
using namespace Geometry;
using namespace Quantities;
using namespace SI;
using namespace UK;
using namespace BIPM;
using namespace Astronomy;
using namespace Constants;
using namespace TestUtilities;

TEST_CLASS(GeometryTest)
{
 public:
  TEST_METHOD(R3ElementTest) {
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
    TestAlternatingBilinearMap(Cross<Speed::Dimensions, Speed::Dimensions>, u,
                               v, w, a, Dimensionless(42));
    TestSymmetricBilinearMap(Dot<Speed::Dimensions, Speed::Dimensions>, u,
                             v, w, a, Dimensionless(42));
  }

  TEST_METHOD(VectorSpaceTests) {
    struct World;
    R3Element<Length> nullDisplacement(0 * Metre, 0 * Metre, 0 * Metre);
    R3Element<Length> u(3 * Metre, -42 * Metre, 0 * Metre);
    R3Element<Length> v(-π * Metre, -e * Metre, -1 * Metre);
    R3Element<Length> w(2 * Metre, 2 * Metre, 2 * Metre);
    R3Element<Length> a(1 * Inch, 2 * Foot, 3 * Admiralty::Fathom);
    TestVectorSpace(Vector<Length, World>(nullDisplacement),
                    Vector<Length, World>(u), Vector<Length, World>(v),
                    Vector<Length, World>(w), Dimensionless(0),
                    Dimensionless(1), Sqrt(163), -Sqrt(2));
    TestVectorSpace(Bivector<Length, World>(nullDisplacement),
                    Bivector<Length, World>(u), Bivector<Length, World>(v),
                    Bivector<Length, World>(w), Dimensionless(0),
                    Dimensionless(1), Sqrt(163), -Sqrt(2));
    TestVectorSpace(Trivector<Length, World>(nullDisplacement.x),
                    Trivector<Length, World>(u.x),
                    Trivector<Length, World>(v.x),
                    Trivector<Length, World>(w.x), Dimensionless(0),
                    Dimensionless(1), Sqrt(163), -Sqrt(2));
  }

};
}

}
