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
    TestAlternatingBilinearMap(Cross<Speed, Speed>, u,
                               v, w, a, Dimensionless(42));
    TestSymmetricPositiveDefiniteBilinearMap(Dot<Speed, Speed>,
                                             u, v, w, a, Dimensionless(42));
  }

  TEST_METHOD(SpecialOrthogonalLieAlgebraTests) {
    struct World;
    R3Element<Dimensionless> u(3, -42, 0);
    R3Element<Dimensionless> v(-π, -e, -1);
    R3Element<Dimensionless> w(2, 2, 2);
    R3Element<Dimensionless> a(1.2, 2.3, 3.4);
    TestLieBracket(Commutator<Dimensionless, Dimensionless, World>,
                   Bivector<Dimensionless, World>(u),
                   Bivector<Dimensionless, World>(v),
                   Bivector<Dimensionless, World>(w),
                   Bivector<Dimensionless, World>(a), Dimensionless(0.42));
  }

  TEST_METHOD(VectorSpaceTests) {
    struct World;
    R3Element<Length> nullDisplacement(0 * Metre, 0 * Metre, 0 * Metre);
    R3Element<Length> u(3 * Metre, -42 * Metre, 0 * Metre);
    R3Element<Length> v(-π * Metre, -e * Metre, -1 * Metre);
    R3Element<Length> w(2 * Metre, 2 * Metre, 2 * Metre);
    R3Element<Length> a(1 * Inch, 2 * Foot, 3 * Admiralty::Fathom);
    {
      std::function<Area(Vector<Length, World>,
                         Vector<Length, World>)> vectorInnerProduct =
      [](Vector<Length, World> a, Vector<Length, World> b) {
        return InnerProduct(a, b);
      };
      std::function<Area(Bivector<Length, World>,
                         Bivector<Length, World>)> bivectorInnerProduct =
      [](Bivector<Length, World> a, Bivector<Length, World> b) {
        return InnerProduct(a, b);
      };
      std::function<Area(Trivector<Length, World>,
                         Trivector<Length, World>)> trivectorInnerProduct =
      [](Trivector<Length, World> a, Trivector<Length, World> b) {
        return InnerProduct(a, b);
      };
      TestInnerProductSpace(vectorInnerProduct,
                            Vector<Length, World>(nullDisplacement),
                            Vector<Length, World>(u), Vector<Length, World>(v),
                            Vector<Length, World>(w), Vector<Length, World>(a),
                            Dimensionless(0), Dimensionless(1), Sqrt(163),
                            -Sqrt(2));
      TestInnerProductSpace(bivectorInnerProduct,
                            Bivector<Length, World>(nullDisplacement),
                            Bivector<Length, World>(u),
                            Bivector<Length, World>(v),
                            Bivector<Length, World>(w),
                            Bivector<Length, World>(a),
                            Dimensionless(0), Dimensionless(1), Sqrt(163),
                            -Sqrt(2));
      TestInnerProductSpace(trivectorInnerProduct,
                            Trivector<Length, World>(nullDisplacement.x),
                            Trivector<Length, World>(u.y),
                            Trivector<Length, World>(v.z),
                            Trivector<Length, World>(w.x),
                            Trivector<Length, World>(a.y),
                            Dimensionless(0), Dimensionless(1), Sqrt(163),
                            -Sqrt(2));
    }
    {
      std::function<Dimensionless(Vector<Dimensionless, World>,
                                  Vector<Dimensionless, World>)>
      vectorInnerProduct = [](Vector<Dimensionless, World> a,
                              Vector<Dimensionless, World> b) {
        return InnerProduct(a, b);
      };
      std::function<Dimensionless(Bivector<Dimensionless, World>,
                                  Bivector<Dimensionless, World>)>
      bivectorInnerProduct = [](Bivector<Dimensionless, World> a,
                                Bivector<Dimensionless, World> b) {
        return InnerProduct(a, b);
      };
      std::function<Dimensionless(Trivector<Dimensionless, World>,
                                  Trivector<Dimensionless, World>)>
      trivectorInnerProduct = [](Trivector<Dimensionless, World> a,
                                 Trivector<Dimensionless, World> b) {
        return InnerProduct(a, b);
      };
      TestInnerProductSpace(vectorInnerProduct,
                            Vector<Dimensionless,
                                   World>(nullDisplacement / Metre),
                            Vector<Dimensionless, World>(u / Metre),
                            Vector<Dimensionless, World>(v / Metre),
                            Vector<Dimensionless, World>(w / Metre),
                            Vector<Dimensionless, World>(a / Metre),
                            Dimensionless(0), Dimensionless(1), Sqrt(163),
                            -Sqrt(2));
      TestInnerProductSpace(bivectorInnerProduct,
                            Bivector<Dimensionless,
                                     World>(nullDisplacement / Metre),
                            Bivector<Dimensionless, World>(u / Metre),
                            Bivector<Dimensionless, World>(v / Metre),
                            Bivector<Dimensionless, World>(w / Metre),
                            Bivector<Dimensionless, World>(a / Metre),
                            Dimensionless(0), Dimensionless(1), Sqrt(163),
                            -Sqrt(2));
      TestInnerProductSpace(trivectorInnerProduct,
                            Trivector<Dimensionless,
                                      World>(nullDisplacement.x / Metre),
                            Trivector<Dimensionless, World>(u.y / Metre),
                            Trivector<Dimensionless, World>(v.z / Metre),
                            Trivector<Dimensionless, World>(w.x / Metre),
                            Trivector<Dimensionless, World>(a.y / Metre),
                            Dimensionless(0), Dimensionless(1), Sqrt(163),
                            -Sqrt(2));
    }
  }
};
}  // namespace GeometryTests
}  // namespace Principia
