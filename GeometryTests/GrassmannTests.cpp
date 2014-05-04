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

TEST_CLASS(GrassmannTests) {
 public:
  struct World;
  TEST_METHOD(SpecialOrthogonalLieAlgebra) {
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

  TEST_METHOD(VectorSpaces) {
    R3Element<Length> nullDisplacement(0 * Metre, 0 * Metre, 0 * Metre);
    R3Element<Length> u(3 * Metre, -42 * Metre, 0 * Metre);
    R3Element<Length> v(-π * Metre, -e * Metre, -1 * Metre);
    R3Element<Length> w(2 * Metre, 2 * Metre, 2 * Metre);
    R3Element<Length> a(1 * Inch, 2 * Foot, 3 * admiralty::Fathom);
    {
      auto vectorInnerProduct = [](Vector<Length, World> a,
                                   Vector<Length, World> b) {
        return InnerProduct(a, b);
      };
      auto bivectorInnerProduct = [](Bivector<Length, World> a,
                                     Bivector<Length, World> b) {
        return InnerProduct(a, b);
      };
      auto trivectorInnerProduct = [](Trivector<Length, World> a,
                                      Trivector<Length, World> b) {
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
      auto vectorInnerProduct = [](Vector<Dimensionless, World> a,
                                   Vector<Dimensionless, World> b) {
        return InnerProduct(a, b);
      };
      auto bivectorInnerProduct = [](Bivector<Dimensionless, World> a,
                                     Bivector<Dimensionless, World> b) {
        return InnerProduct(a, b);
      };
      auto trivectorInnerProduct = [](Trivector<Dimensionless, World> a,
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

  TEST_METHOD(GrassmannAlgebra) {
    R3Element<Dimensionless> u(3, -42, 0);
    R3Element<Dimensionless> v(-π, -e, -1);
    R3Element<Dimensionless> w(2, 2, 2);
    R3Element<Dimensionless> a(1.2, 2.3, 3.4);
    auto vectorWedge = [](Vector<Dimensionless, World> a,
                          Vector<Dimensionless, World> b) {
        return Wedge(a, b);
      };
    TestAlternatingBilinearMap(vectorWedge, Vector<Dimensionless, World>(u),
                               Vector<Dimensionless, World>(u),
                               Vector<Dimensionless, World>(u),
                               Vector<Dimensionless, World>(u),
                               Dimensionless(6 * 9));
  }
};

}  // namespace geometry
}  // namespace principia
