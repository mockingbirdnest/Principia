#include "geometry/grassmann.hpp"

#include <cfloat>
#include <functional>
#include <iostream>  // NOLINT(readability/streams)
#include <utility>

#include "base/cpuid.hpp"
#include "geometry/frame.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/fma.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/constants.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/algebra.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/check_well_formedness.hpp"  // 🧙 For PRINCIPIA_CHECK_ILL_FORMED.

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::base::_cpuid;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::quantities::_uk;
using namespace principia::quantities::_uk::admiralty;
using namespace principia::testing_utilities::_algebra;
using namespace principia::testing_utilities::_almost_equals;

struct TransparentInnerProduct final {
  template<typename Left, typename Right>
  decltype(InnerProduct(std::forward<Left>(std::declval<Left&&>()),
                        std::forward<Right>(std::declval<Right&&>())))
  operator()(Left&& left, Right&& right)  // NOLINT(whitespace/operators)
      const {
    return InnerProduct(std::forward<Left>(left),
                        std::forward<Right>(right));
  }
};

struct TransparentWedge final {
  template<typename Left, typename Right>
  decltype(Wedge(std::forward<Left>(std::declval<Left&&>()),
                 std::forward<Right>(std::declval<Right&&>())))
  operator()(Left&& left, Right&& right)  // NOLINT(whitespace/operators)
      const {
    return Wedge(std::forward<Left>(left),
                 std::forward<Right>(right));
  }
};

class GrassmannTest : public testing::Test {
 public:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

 protected:
  R3Element<Length> const null_displacement_ = {0 * Metre,
                                                0 * Metre,
                                                0 * Metre};
  R3Element<Length> const u_ = {3 * Metre, -42 * Metre, 0 * Metre};
  R3Element<Length> const v_ = {-π * Metre, -e * Metre, -1 * Metre};
  R3Element<Length> const w_ = {2 * Metre, 2 * Metre, 2 * Metre};
  R3Element<Length> const a_ = {1 * Inch, 2 * Foot, 3 * Fathom};
};

using GrassmannDeathTest = GrassmannTest;

TEST_F(GrassmannTest, VectorFMA) {
  if (!CanEmitFMAInstructions || !CPUIDFeatureFlag::FMA.IsSet()) {
    GTEST_SKIP() << "Cannot test FMA on a machine without FMA";
  }
  Length const a = a_.x;
  Vector<Length, World> const v(v_);
  Vector<Area, World> const w(a_.y * w_);
  EXPECT_THAT(FusedMultiplyAdd(a, v, w), AlmostEquals(a * v + w, 0));
  EXPECT_THAT(FusedMultiplySubtract(a, v, w), AlmostEquals(a * v - w, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(a, v, w), AlmostEquals(-a * v + w, 0));
  EXPECT_THAT(FusedNegatedMultiplySubtract(a, v, w),
              AlmostEquals(-a * v - w, 0));
  EXPECT_THAT(FusedMultiplyAdd(v, a, w), AlmostEquals(v * a + w, 0));
  EXPECT_THAT(FusedMultiplySubtract(v, a, w), AlmostEquals(v * a - w, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(v, a, w), AlmostEquals(-v * a + w, 0));
  EXPECT_THAT(FusedNegatedMultiplySubtract(v, a, w),
              AlmostEquals(-v * a - w, 0));
}

TEST_F(GrassmannTest, BivectorFMA) {
  if (!CanEmitFMAInstructions || !CPUIDFeatureFlag::FMA.IsSet()) {
    GTEST_SKIP() << "Cannot test FMA on a machine without FMA";
  }
  Length const a = a_.x;
  Bivector<Length, World> const v(v_);
  Bivector<Area, World> const w(a_.y * w_);
  EXPECT_THAT(FusedMultiplyAdd(a, v, w), AlmostEquals(a * v + w, 0));
  EXPECT_THAT(FusedMultiplySubtract(a, v, w), AlmostEquals(a * v - w, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(a, v, w), AlmostEquals(-a * v + w, 0));
  EXPECT_THAT(FusedNegatedMultiplySubtract(a, v, w),
              AlmostEquals(-a * v - w, 0));
  EXPECT_THAT(FusedMultiplyAdd(v, a, w), AlmostEquals(v * a + w, 0));
  EXPECT_THAT(FusedMultiplySubtract(v, a, w), AlmostEquals(v * a - w, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(v, a, w), AlmostEquals(-v * a + w, 0));
  EXPECT_THAT(FusedNegatedMultiplySubtract(v, a, w),
              AlmostEquals(-v * a - w, 0));
}

TEST_F(GrassmannTest, TrivectorFMA) {
  if (!CanEmitFMAInstructions || !CPUIDFeatureFlag::FMA.IsSet()) {
    GTEST_SKIP() << "Cannot test FMA on a machine without FMA";
  }
  Length const a = a_.x;
  Trivector<Length, World> const v(v_.x);
  Trivector<Area, World> const w(a_.y * w_.x);
  EXPECT_THAT(FusedMultiplyAdd(a, v, w), AlmostEquals(a * v + w, 0));
  EXPECT_THAT(FusedMultiplySubtract(a, v, w), AlmostEquals(a * v - w, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(a, v, w), AlmostEquals(-a * v + w, 0));
  EXPECT_THAT(FusedNegatedMultiplySubtract(a, v, w),
              AlmostEquals(-a * v - w, 0));
  EXPECT_THAT(FusedMultiplyAdd(v, a, w), AlmostEquals(v * a + w, 0));
  EXPECT_THAT(FusedMultiplySubtract(v, a, w), AlmostEquals(v * a - w, 0));
  EXPECT_THAT(FusedNegatedMultiplyAdd(v, a, w), AlmostEquals(-v * a + w, 0));
  EXPECT_THAT(FusedNegatedMultiplySubtract(v, a, w),
              AlmostEquals(-v * a - w, 0));
}

TEST_F(GrassmannTest, Operators) {
  TestEquality(Bivector<Length, World>(u_),
                                  Bivector<Length, World>(v_));
  TestEquality(Vector<Length, World>(u_),
                                  Vector<Length, World>(v_));
  TestEquality(Trivector<Length, World>(u_.x),
                                  Trivector<Length, World>(v_.x));
}

TEST_F(GrassmannTest, SpecialOrthogonalLieAlgebra) {
  TestLieBracket(Commutator<double, double, World>,
                 Bivector<double, World>(u_ / Foot),
                 Bivector<double, World>(v_ / Metre),
                 Bivector<double, World>(w_ / Rod),
                 Bivector<double, World>(a_ / Furlong),
                 0.42, 0, 1);
}

TEST_F(GrassmannTest, MixedScalarMultiplication) {
  TestBilinearMap(std::multiplies<>(),
                  1 / Second,
                  1 / JulianYear,
                  Vector<Length, World>(u_),
                  Vector<Length, World>(v_),
                  42, 0, 1);
  TestBilinearMap(std::multiplies<>(),
                  Vector<Length, World>(w_),
                  Vector<Length, World>(a_),
                  -1 / Day,
                  SpeedOfLight / Parsec,
                  -π, 0, 1);
  Inverse<Time> t = -3 / Second;
  EXPECT_EQ((t * Vector<Length, World>(u_)), (Vector<Length, World>(u_) * t));
  EXPECT_EQ((Vector<Length, World>(v_) * t) / t, (Vector<Length, World>(v_)));
}

TEST_F(GrassmannTest, VectorSpaces) {
  TestInnerProductSpace(TransparentInnerProduct(),
                        Vector<Length, World>(null_displacement_),
                        Vector<Length, World>(u_),
                        Vector<Length, World>(v_),
                        Vector<Length, World>(w_),
                        Vector<Length, World>(a_),
                        0.0, 1.0, Sqrt(163), -Sqrt(2), 0, 18);
  TestInnerProductSpace(TransparentInnerProduct(),
                        Bivector<Length, World>(null_displacement_),
                        Bivector<Length, World>(u_),
                        Bivector<Length, World>(v_),
                        Bivector<Length, World>(w_),
                        Bivector<Length, World>(a_),
                        0.0, 1.0, Sqrt(163), -Sqrt(2), 0,  18);
  TestInnerProductSpace(TransparentInnerProduct(),
                        Trivector<Length, World>(null_displacement_.x),
                        Trivector<Length, World>(u_.y),
                        Trivector<Length, World>(v_.z),
                        Trivector<Length, World>(w_.x),
                        Trivector<Length, World>(a_.y),
                        0.0, 1.0, Sqrt(163), -Sqrt(2), 0, 1);
  TestInnerProductSpace(TransparentInnerProduct(),
                        Vector<double, World>(null_displacement_ / Metre),
                        Vector<double, World>(u_ / Metre),
                        Vector<double, World>(v_ / Metre),
                        Vector<double, World>(w_ / Metre),
                        Vector<double, World>(a_ / Metre),
                        0.0, 1.0, Sqrt(163), -Sqrt(2), 0, 18);
  TestInnerProductSpace(TransparentInnerProduct(),
                        Bivector<double, World>(null_displacement_ / Metre),
                        Bivector<double, World>(u_ / Metre),
                        Bivector<double, World>(v_ / Metre),
                        Bivector<double, World>(w_ / Metre),
                        Bivector<double, World>(a_ / Metre),
                        0.0, 1.0, Sqrt(163), -Sqrt(2), 0, 18);
  TestInnerProductSpace(TransparentInnerProduct(),
                        Trivector<double, World>(null_displacement_.x / Metre),
                        Trivector<double, World>(u_.y / Metre),
                        Trivector<double, World>(v_.z / Metre),
                        Trivector<double, World>(w_.x / Metre),
                        Trivector<double, World>(a_.y / Metre),
                        0.0, 1.0, Sqrt(163), -Sqrt(2), 0, 1);
}

TEST_F(GrassmannTest, GrassmannAlgebra) {
  TestAlternatingBilinearMap(TransparentWedge(),
                             Vector<double, World>(u_ / Metre),
                             Vector<double, World>(v_ / Metre),
                             Vector<double, World>(w_ / Metre),
                             Vector<double, World>(a_ / Metre),
                             6.0 * 9.0, 0, 1);
  TestBilinearMap(TransparentWedge(),
                  Vector<Length, World>(u_),
                  Vector<Length, World>(v_),
                  Bivector<Speed, World>(w_ / Second),
                  Bivector<Speed, World>(a_ / Second),
                  6.0 * 9.0, 0, 1);
  TestBilinearMap(TransparentWedge(),
                  Bivector<Length, World>(u_),
                  Bivector<Length, World>(v_),
                  Vector<Speed, World>(w_ / Second),
                  Vector<Speed, World>(a_ / Second),
                  6.0 * 9.0, 0, 1);
  EXPECT_EQ(
      Wedge(Vector<Speed, World>(v_ / Second), Bivector<Length, World>(u_)),
      Wedge(Bivector<Length, World>(u_), Vector<Speed, World>(v_ / Second)));
}

// The Greek letters cause a warning when stringified by the macros, because
// apparently Visual Studio doesn't encode strings in UTF-8 by default.
#pragma warning(disable: 4566)

TEST_F(GrassmannTest, Actions) {
  Vector<Length, World> const a(u_);
  Vector<Length, World> const b(v_);
  Bivector<Length, World> const β(v_);
  Bivector<Length, World> const γ(w_);
  // A strongly typed version of the Lagrange formula
  // a × (b × c) = b (a · c) - c (a · b).
  EXPECT_THAT(a * Commutator(β, γ),
              AlmostEquals(β * Wedge(a, γ) - γ * Wedge(a, β), 26));
  EXPECT_THAT(Commutator(β, γ) * a,
              AlmostEquals(Wedge(a, β) * γ - β * Wedge(a, γ), 26));

  EXPECT_THAT(a * Wedge(b, γ), AlmostEquals(Wedge(γ, b) * a, 0, 21));
}

#pragma warning(default: 4566)

TEST_F(GrassmannTest, Norm) {
  Vector<Length, World> const v({-3 * 4 * Metre, 4 * 4 * Metre, 5 * 3 * Metre});
  EXPECT_THAT(v.Norm(), Eq(5 * 5 * Metre));
  Bivector<Length, World> const w({+20 * 21 * Metre,
                                   -21 * 21 * Metre,
                                   +29 * 20 * Metre});
  EXPECT_THAT(w.Norm(), Eq(29 * 29 * Metre));
  Trivector<Length, World> const u(-4 * Furlong);
  EXPECT_THAT(u.Norm(), Eq(4 * Furlong));
}

TEST_F(GrassmannTest, Normalize) {
  Vector<Length, World> const v({-1 * Metre, 2 * Metre, 3 * Metre});
  EXPECT_THAT(Normalize(v),
              Eq(Vector<double, World>({-1 / Sqrt(14),
                                        2 / Sqrt(14),
                                        3 / Sqrt(14)})));
  Bivector<Length, World> const w({4 * Metre, -5 * Metre, 6 * Metre});
  EXPECT_THAT(Normalize(w),
              Eq(Bivector<double, World>({4 / Sqrt(77),
                                          -5 / Sqrt(77),
                                          6 / Sqrt(77)})));
  Trivector<Length, World> const u(-4 * Furlong);
  EXPECT_THAT(Normalize(u), Eq(Trivector<double, World>(-1)));
}

// Check that non-serializable frames are detected at compile-time.

using F = Frame<struct FrameTag>;

PRINCIPIA_CHECK_WELL_FORMED(
    v.WriteToMessage(&message),
    with_variable<Vector<Length, GrassmannTest::World>> v,
    with_variable<serialization::Multivector> message);
// TODO(phl): We should refuse to serialize these at compile time; right now
// only deserialization fails.
PRINCIPIA_CHECK_WELL_FORMED(v.WriteToMessage(&message),
                            with_variable<Vector<Length, F>> v,
                            with_variable<serialization::Multivector> message);
PRINCIPIA_CHECK_WELL_FORMED(v.WriteToMessage(&message),
                            with_variable<Bivector<Length, F>> v,
                            with_variable<serialization::Multivector> message);
PRINCIPIA_CHECK_WELL_FORMED(v.WriteToMessage(&message),
                            with_variable<Trivector<Length, F>> v,
                            with_variable<serialization::Multivector> message);

PRINCIPIA_CHECK_WELL_FORMED_WITH_TYPES(
    V::ReadFromMessage(message),
    (typename V = Vector<Length, GrassmannTest::World>),
    with_variable<serialization::Multivector> message);
PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(
    V::ReadFromMessage(message),
    (typename V = Vector<Length, F>),
    with_variable<serialization::Multivector> message);
PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(
    V::ReadFromMessage(message),
    (typename V = Bivector<Length, F>),
    with_variable<serialization::Multivector> message);
PRINCIPIA_CHECK_ILL_FORMED_WITH_TYPES(
    V::ReadFromMessage(message),
    (typename V = Trivector<Length, F>),
    with_variable<serialization::Multivector> message);

TEST_F(GrassmannDeathTest, SerializationError) {
  using V = Vector<Length, World>;
  using B = Bivector<Length, World>;
  using T = Trivector<Length, World>;

  EXPECT_DEATH({
    serialization::Multivector message;
    V const v({-1 * Metre, 2 * Metre, 3 * Metre});
    v.WriteToMessage(&message);
    [[maybe_unused]] B const b = B::ReadFromMessage(message);
  }, "has_bivector");
  EXPECT_DEATH({
    serialization::Multivector message;
    B const b({-1 * Metre, 2 * Metre, 3 * Metre});
    b.WriteToMessage(&message);
    [[maybe_unused]] T const t = T::ReadFromMessage(message);
  }, "has_trivector");
  EXPECT_DEATH({
    serialization::Multivector message;
    T const t(1 * Metre);
    t.WriteToMessage(&message);
    [[maybe_unused]] V const v = V::ReadFromMessage(message);
  }, "has_vector");
}

TEST_F(GrassmannTest, SerializationSuccess) {
  serialization::Multivector message;

  Vector<Length, World> const v({-1 * Metre, 2 * Metre, 3 * Metre});
  v.WriteToMessage(&message);
  EXPECT_TRUE(message.has_vector());
  EXPECT_FALSE(message.has_bivector());
  EXPECT_FALSE(message.has_trivector());
  Vector<Length, World> const w =
      Vector<Length, World>::ReadFromMessage(message);
  EXPECT_EQ(v, w);

  Bivector<Pressure, World> const b({-4 * Pascal, -5 * Pascal, 6 * Pascal});
  b.WriteToMessage(&message);
  EXPECT_FALSE(message.has_vector());
  EXPECT_TRUE(message.has_bivector());
  EXPECT_FALSE(message.has_trivector());
  Bivector<Pressure, World> const c =
      Bivector<Pressure, World>::ReadFromMessage(message);
  EXPECT_EQ(b, c);

  Trivector<Charge, World> const t(-7 * Coulomb);
  t.WriteToMessage(&message);
  EXPECT_FALSE(message.has_vector());
  EXPECT_FALSE(message.has_bivector());
  EXPECT_TRUE(message.has_trivector());
  EXPECT_EQ(-7, message.trivector().magnitude());
  Trivector<Charge, World> const u =
      Trivector<Charge, World>::ReadFromMessage(message);
  EXPECT_EQ(t, u);
}

TEST_F(GrassmannTest, Angles) {
  Vector<double, World> const x({1, 0, 0});
  Vector<double, World> const y({0, 1, 0});
  Vector<double, World> const z({0, 0, 1});
  EXPECT_THAT(AngleBetween(x, y), Eq(π / 2 * Radian));
  EXPECT_THAT(AngleBetween(x, -y), Eq(π / 2 * Radian));
  EXPECT_THAT(AngleBetween(x, x + y), Eq(π / 4 * Radian));
  EXPECT_THAT(OrientedAngleBetween(x, y, Wedge(x, y)), Eq(π / 2 * Radian));
  EXPECT_THAT(OrientedAngleBetween(x, -y, Wedge(x, y)), Eq(-π / 2 * Radian));
  EXPECT_THAT(AngleBetween(Wedge(x, y), Wedge(y, z)), Eq(π / 2 * Radian));
}

}  // namespace geometry
}  // namespace principia
