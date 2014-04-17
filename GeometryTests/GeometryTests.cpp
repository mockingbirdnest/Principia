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
    TestBilinearMap(Cross<Speed::Dimensions, Speed::Dimensions>, u, v, w, a,
                    Dimensionless(42));
  }

  TEST_METHOD(VectorSpaceTests) {
    struct World;
    Vector<Length, World> nullWorldDisplacement(R3Element<Length>(0 * Metre,
                                                                  0 * Metre,
                                                                  0 * Metre));
    Vector<Length, World> u(R3Element<Length>(3 * Metre, 
                                              -42 * Metre,
                                              0 * Metre));
    Vector<Length, World> v(R3Element<Length>(-π * Metre,
                                              -e * Metre,
                                              -1 * Metre));
    Vector<Length, World> w(R3Element<Length>(2 * Metre,
                                              2 * Metre,
                                              2 * Metre));
    TestVectorSpace(nullWorldDisplacement, u, v, w, Dimensionless(0),
                    Dimensionless(1), Sqrt(163), -Sqrt(2));
  }

};
}

}
