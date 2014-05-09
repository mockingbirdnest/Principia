#include "stdafx.hpp"

#include <float.h>

#include <CppUnitTest.h>
#include <functional>

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
#include "TestUtilities/ExplicitOperators.hpp"
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
  R3Element<Length> const null_displacement_ =
      R3Element<Length>(0 * Metre, 0 * Metre, 0 * Metre);
  R3Element<Length> const u_ =
      R3Element<Length>(3 * Metre, -42 * Metre, 0 * Metre);
  R3Element<Length> const v_ =
      R3Element<Length>(-π * Metre, -e * Metre, -1 * Metre);
  R3Element<Length> const w_ =
      R3Element<Length>(2 * Metre, 2 * Metre, 2 * Metre);
  R3Element<Length> const a_ =
      R3Element<Length>(1 * Inch, 2 * Foot, 3 * admiralty::Fathom);

  template<typename LScalar, typename RScalar, typename Frame, int Rank>
  static Product<LScalar, RScalar> MultivectorInnerProduct(
      Multivector<LScalar, Frame, Rank> const& left,
      Multivector<RScalar, Frame, Rank> const& right) {
    return InnerProduct(left, right);
  }

  template<typename LScalar, typename RScalar, typename Frame, int LRank,
          int RRank>
  static Multivector<Product<LScalar, RScalar>, Frame, LRank + RRank>
  Multivectorwedge(
      Multivector<LScalar, Frame, LRank> const& left,
      Multivector<RScalar, Frame, RRank> const& right) {
    return Wedge(left, right);
  }

  TEST_METHOD(SpecialOrthogonalLieAlgebra) {
    TestLieBracket(Commutator<Dimensionless, Dimensionless, World>,
                   Bivector<Dimensionless, World>(u_ / Foot),
                   Bivector<Dimensionless, World>(v_ / Metre),
                   Bivector<Dimensionless, World>(w_ / Rod),
                   Bivector<Dimensionless, World>(a_ / Furlong),
                   Dimensionless(0.42), 128 * DBL_EPSILON);
  }

  TEST_METHOD(MixedScalarMultiplication) {
     TestBilinearMap(
        Times<Vector<Speed, World>, Inverse<Time>, Vector<Length, World>>,
        1 / Second,
        1 / JulianYear,
        Vector<Length, World>(u_),
        Vector<Length, World>(v_),
        Dimensionless(42),
        2 * DBL_EPSILON);
     TestBilinearMap(
        Times<Vector<Speed, World>, Vector<Length, World>, Inverse<Time>>,
        Vector<Length, World>(w_),
        Vector<Length, World>(a_),
        -1 / Day,
        SpeedOfLight / Parsec,
        Dimensionless(-π),
        DBL_EPSILON);
     Inverse<Time> t = -3 / Second;
     AssertEqual(t * Vector<Length, World>(u_), Vector<Length, World>(u_) * t);
     AssertEqual((Vector<Length, World>(v_) * t) / t,
                 Vector<Length, World>(v_));
  }

  TEST_METHOD(VectorSpaces) {
    TestInnerProductSpace(MultivectorInnerProduct<Length, Length, World, 1>,
                          Vector<Length, World>(null_displacement_),
                          Vector<Length, World>(u_),
                          Vector<Length, World>(v_),
                          Vector<Length, World>(w_),
                          Vector<Length, World>(a_),
                          Dimensionless(0), Dimensionless(1), Sqrt(163),
                          -Sqrt(2), 12 * DBL_EPSILON);
    TestInnerProductSpace(MultivectorInnerProduct<Length, Length, World, 2>,
                          Bivector<Length, World>(null_displacement_),
                          Bivector<Length, World>(u_),
                          Bivector<Length, World>(v_),
                          Bivector<Length, World>(w_),
                          Bivector<Length, World>(a_),
                          Dimensionless(0), Dimensionless(1), Sqrt(163),
                          -Sqrt(2), 12 * DBL_EPSILON);
    TestInnerProductSpace(MultivectorInnerProduct<Length, Length, World, 3>,
                          Trivector<Length, World>(null_displacement_.x),
                          Trivector<Length, World>(u_.y),
                          Trivector<Length, World>(v_.z),
                          Trivector<Length, World>(w_.x),
                          Trivector<Length, World>(a_.y),
                          Dimensionless(0), Dimensionless(1), Sqrt(163),
                          -Sqrt(2));
    TestInnerProductSpace(MultivectorInnerProduct<Dimensionless, Dimensionless,
                                                  World, 1>,
                          Vector<Dimensionless,
                                 World>(null_displacement_ / Metre),
                          Vector<Dimensionless, World>(u_ / Metre),
                          Vector<Dimensionless, World>(v_ / Metre),
                          Vector<Dimensionless, World>(w_ / Metre),
                          Vector<Dimensionless, World>(a_ / Metre),
                          Dimensionless(0), Dimensionless(1), Sqrt(163),
                          -Sqrt(2), 12 * DBL_EPSILON);
    TestInnerProductSpace(MultivectorInnerProduct<Dimensionless, Dimensionless,
                                                  World, 2>,
                          Bivector<Dimensionless,
                                   World>(null_displacement_ / Metre),
                          Bivector<Dimensionless, World>(u_ / Metre),
                          Bivector<Dimensionless, World>(v_ / Metre),
                          Bivector<Dimensionless, World>(w_ / Metre),
                          Bivector<Dimensionless, World>(a_ / Metre),
                          Dimensionless(0), Dimensionless(1), Sqrt(163),
                          -Sqrt(2), 12 * DBL_EPSILON);
    TestInnerProductSpace(MultivectorInnerProduct<Dimensionless, Dimensionless,
                                                  World, 3>,
                          Trivector<Dimensionless,
                                    World>(null_displacement_.x / Metre),
                          Trivector<Dimensionless, World>(u_.y / Metre),
                          Trivector<Dimensionless, World>(v_.z / Metre),
                          Trivector<Dimensionless, World>(w_.x / Metre),
                          Trivector<Dimensionless, World>(a_.y / Metre),
                          Dimensionless(0), Dimensionless(1), Sqrt(163),
                          -Sqrt(2));
  }

  TEST_METHOD(GrassmannAlgebra) {
    TestAlternatingBilinearMap(
        Multivectorwedge<Dimensionless, Dimensionless, World, 1, 1>,
        Vector<Dimensionless, World>(u_ / Metre),
        Vector<Dimensionless, World>(v_ / Metre),
        Vector<Dimensionless, World>(w_ / Metre),
        Vector<Dimensionless, World>(a_ / Metre),
        Dimensionless(6 * 9),
        2 * DBL_EPSILON);
  }
};

}  // namespace geometry
}  // namespace principia
