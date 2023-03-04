#include "numerics/double_precision.hpp"

#include <limits>
#include <random>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

#define PRINCIPIA_USE_IACA 0
#if PRINCIPIA_USE_IACA
#include "intel/iacaMarks.h"
#endif

namespace principia {
namespace numerics {

using testing_utilities::AlmostEquals;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Ne;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_named_quantities;
using namespace principia::geometry::_point;
using namespace principia::geometry::_r3_element;
using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

constexpr double ε = std::numeric_limits<double>::epsilon();
constexpr double ε² = ε * ε;
constexpr double ε³ = ε² * ε;
constexpr double ε⁴ = ε³ * ε;

using World = Frame<struct WorldTag>;

class DoublePrecisionTest : public ::testing::Test {};

#if PRINCIPIA_USE_IACA
// A convenient skeleton for analysing code with IACA.
TEST_F(DoublePrecisionTest, DISABLED_IACA) {
  auto iaca = [](DoublePrecision<Displacement<World>> const& x,
                 DoublePrecision<Displacement<World>> const& y) {
    IACA_VC64_START;
    auto const result = x - y;
    IACA_VC64_END;
    return result;
  };
  auto const x = Displacement<World>({1 * Metre, 2 * Metre, 3 * Metre});
  auto const y = Displacement<World>({4 * Metre, 3 * Metre, 2 * Metre});
  CHECK_EQ(iaca(DoublePrecision<Displacement<World>>(x),
                DoublePrecision<Displacement<World>>(y)),
           iaca(DoublePrecision<Displacement<World>>(x),
                DoublePrecision<Displacement<World>>(y)));
}
#endif

TEST_F(DoublePrecisionTest, CompensatedSummationIncrement) {
  Position<World> const initial =
      World::origin + Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre});
  DoublePrecision<Position<World>> accumulator(initial);
  Displacement<World> const δ({ε / 4 * Metre, 0 * Metre, 0 * Metre});
  for (int i = 0; i < 4; ++i) {
    accumulator.Increment(δ);
    if (i < 2) {
      EXPECT_THAT(accumulator.value, Eq(initial));
    }
  }
  EXPECT_THAT((accumulator.value - World::origin).coordinates().x,
              Eq((1 + ε) * Metre));
  EXPECT_THAT(accumulator.error.coordinates().x, Eq(0 * Metre));
}

TEST_F(DoublePrecisionTest, CompensatedSummationDecrement) {
  Position<World> const initial =
      World::origin + Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre});
  DoublePrecision<Position<World>> accumulator(initial);
  Displacement<World> const δ({ε / 4 * Metre, 0 * Metre, 0 * Metre});
  for (int i = 0; i < 4; ++i) {
    accumulator.Decrement(δ);
    if (i < 1) {
      EXPECT_THAT(accumulator.value, Eq(initial));
    }
  }
  EXPECT_THAT((accumulator.value - World::origin).coordinates().x,
              Eq((1 - ε) * Metre));
  EXPECT_THAT(accumulator.error.coordinates().x, Eq(0 * Metre));
}

TEST_F(DoublePrecisionTest, IllConditionedCompensatedSummationIncrement) {
  Length const x = (1 + ε) * Metre;
  Point<Length> const zero;
  for (bool const cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Point<Length>> accumulator;
    accumulator.Increment(+x);
    accumulator.Increment(-y);
    accumulator.Increment(+x);
    accumulator.Increment(-y);
    accumulator.Increment(-x);
    accumulator.Increment(+y);
    accumulator.Increment(-x);
    accumulator.Increment(+y);
    if (cancellation) {
      EXPECT_THAT(accumulator.value - zero, Eq(ε² * Metre));
    } else {
      EXPECT_THAT(accumulator.value - zero, Eq(0 * Metre));
    }
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

TEST_F(DoublePrecisionTest, IllConditionedCompensatedSummationDecrement) {
  Length const x = (1 + ε) * Metre;
  Point<Length> const zero;
  for (bool const cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Point<Length>> accumulator;
    accumulator.Decrement(+x);
    accumulator.Decrement(-y);
    accumulator.Decrement(+x);
    accumulator.Decrement(-y);
    accumulator.Decrement(-x);
    accumulator.Decrement(+y);
    accumulator.Decrement(-x);
    accumulator.Decrement(+y);
    if (cancellation) {
      EXPECT_THAT(accumulator.value - zero, Eq(-ε² * Metre));
    } else {
      EXPECT_THAT(accumulator.value - zero, Eq(0 * Metre));
    }
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

TEST_F(DoublePrecisionTest, LongAdd) {
  Length const x = (1 + ε) * Metre;
  Point<Length> const zero;
  for (bool cancellation : {true, false}) {
    Length const y = cancellation ? ε * x : π * x;
    DoublePrecision<Point<Length>> accumulator;
    accumulator += TwoSum(+x, -y);
    accumulator -= TwoSum(-x, +y);
    accumulator -= TwoSum(+x, -y);
    accumulator += TwoSum(-x, +y);
    EXPECT_THAT(accumulator.value - zero, Eq(0 * Metre));
    EXPECT_THAT(accumulator.error, Eq(0 * Metre));
  }
}

TEST_F(DoublePrecisionTest, LongAddPositions) {
  DoublePrecision<Position<World>> accumulator;
  Displacement<World> const δ_value(
      {1 * Metre, 2 * Metre, 3 * Metre});
  Displacement<World> const δ_error(
      {ε / 4 * Metre, ε / 2 * Metre, ε / 2 * Metre});
  DoublePrecision<Displacement<World>> const δ = TwoSum(δ_error, δ_value);
  EXPECT_THAT(δ.value, Eq(δ_value));
  EXPECT_THAT(δ.error, Eq(δ_error));
  for (int i = 0; i < 4; ++i) {
    accumulator += δ;
  }
  for (int i = 0; i < 3; ++i) {
    accumulator -= DoublePrecision<Displacement<World>>(δ_value);
  }
  DoublePrecision<Displacement<World>> const accumulated_displacement =
      accumulator - DoublePrecision<Position<World>>(World::origin);
  EXPECT_THAT(accumulated_displacement.value,
              Eq(δ_value + Displacement<World>(
                               {ε * Metre, 2 * ε * Metre, 2 * ε * Metre})));
  EXPECT_THAT(accumulated_displacement.error,
              Eq(Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre})));
}

TEST_F(DoublePrecisionTest, DoubleDoubleDouble) {
  using DoubleDouble = DoublePrecision<double>;
  using DoubleDoubleDouble = DoublePrecision<DoubleDouble>;
  auto wider = [](double x) { return DoubleDouble(x); };
  auto widest = [](double x) { return DoubleDoubleDouble(DoubleDouble(x)); };
  DoubleDoubleDouble accumulator;
  accumulator.Increment(wider(1));
  accumulator.Increment(wider(ε / 2));
  DoubleDoubleDouble const δ(wider(ε² / 4));
  EXPECT_THAT(DebugString(accumulator + δ),
              Eq("+1.00000000000000000e+00|+1.11022302462515654e-16|"
                 "+1.23259516440783095e-32|+0.00000000000000000e+00"));
  for (int i = 0; i < 4; ++i) {
    accumulator += δ;
  }
  accumulator -= widest(ε / 2);
  EXPECT_THAT(DebugString(accumulator),
              Eq("+1.00000000000000000e+00|+4.93038065763132378e-32|"
                 "+0.00000000000000000e+00|+0.00000000000000000e+00"));
  DoubleDoubleDouble const long_long_long_float =
      TwoSum(TwoSum(1, ε / 2), TwoSum(ε² / 4, ε³ / 8));
  EXPECT_THAT(DebugString(long_long_long_float),
              Eq("+1.00000000000000000e+00|+1.11022302462515654e-16|"
                 "+1.23259516440783095e-32|+1.36845553156720417e-48"));
  EXPECT_THAT(long_long_long_float + widest(ε⁴ / 16), Eq(long_long_long_float));
  EXPECT_THAT(long_long_long_float + long_long_long_float,
              Eq(TwoSum(Scale(2, long_long_long_float.value),
                        Scale(2, long_long_long_float.error))));
}

TEST_F(DoublePrecisionTest, ComparableTwoSum) {
  DoublePrecision<double> x = TwoSum(π, e);
  int value_exponent;
  int error_exponent;
  std::frexp(x.value, &value_exponent);
  std::frexp(x.error, &error_exponent);
  EXPECT_THAT(value_exponent - error_exponent, Ge(53));
}

// There is a some of replicated code (up to signs) in the differences because
// of typing concerns.  We check consistency of many things with the operator+
// on DoublePrecision<Vector>.
TEST_F(DoublePrecisionTest, Consistencies) {
  using Vector = R3Element<double>;
  using Point = Point<Vector>;
  DoublePrecision<Point> const wide_origin{};
  Vector const null_vector{0, 0, 0};
  Vector const v1{π, -e, Sqrt(2)};
  Vector const v2 = 1024 * ε * Vector{(1 + Sqrt(5)) / 2, std::log(2), -Sqrt(π)};
  Vector const v3{std::log(π), Tan(e * Radian), Sqrt(e)};
  Vector const v4 =
      1024 * ε *
      Vector{std::log((1 + Sqrt(5)) / 2), Sqrt(std::log(2)), std::log(Sqrt(2))};
  DoublePrecision<Vector> const wide_v1(v1);
  DoublePrecision<Vector> const w1 = TwoSum(v1, v2);
  DoublePrecision<Vector> const w2 = TwoSum(v3, v4);
  DoublePrecision<Point> const q1 = wide_origin + w1;
  DoublePrecision<Point> const q2 = wide_origin + w2;

  // DoublePrecision<Vector> - DoublePrecision<Vector>.
  EXPECT_THAT(DebugString(w1 - w2), Eq(DebugString(w1 + -w2)));
  // DoublePrecision<Point> - DoublePrecision<Vector>.
  EXPECT_THAT(DebugString((q1 - w2) - wide_origin), Eq(DebugString(w1 + -w2)));
  // DoublePrecision<Point> - DoublePrecision<Point>.
  EXPECT_THAT(DebugString(q1 - q2), Eq(DebugString(w1 + -w2)));

  // We now showcase the difference between |Increment| and |operator+=|.
  auto compensated_accumulator = -w2;
  compensated_accumulator.Increment(v1);
  auto double_accumulator = -w2;
  double_accumulator += wide_v1;
  EXPECT_THAT(compensated_accumulator.value,
              AlmostEquals(double_accumulator.value, 1));
  EXPECT_THAT(compensated_accumulator.error, Eq(null_vector));
  EXPECT_THAT(double_accumulator.error, Ne(null_vector));
  compensated_accumulator.Increment(-v1);
  double_accumulator -= wide_v1;
  EXPECT_THAT(double_accumulator, Eq(-w2));
  EXPECT_THAT(compensated_accumulator, Ne(-w2));
}

TEST_F(DoublePrecisionTest, Product) {
  Mass const a = 1.0 / 3.0 * Kilogram;
  Speed const b = 1.0 / 7.0 * Metre / Second;
  DoublePrecision<Momentum> const c = TwoProduct(a, b);
  DoublePrecision<Momentum> const d = VeltkampDekkerProduct(a, b);
  EXPECT_THAT(c.value, AlmostEquals(d.value, 0));
  EXPECT_THAT(c.error, AlmostEquals(d.error, 0));
  // The numbers below were obtained using Mathematica.
  EXPECT_THAT(c.value,
              AlmostEquals(6862628003612184.0 * 0x1p-57 * Kilogram *
                               Metre / Second,
                           0));
  EXPECT_THAT(c.error,
              AlmostEquals(-3431314001806092.0 * 0x1p-110 * Kilogram *
                               Metre / Second,
                           0));
}

TEST_F(DoublePrecisionTest, LongProduct) {
  DoublePrecision<Length> a(3 * Metre);
  a.Increment(474 * ε * Metre);
  DoublePrecision<Length> b(7 * Metre);
  b.Increment(-517 * ε * Metre);
  DoublePrecision<Square<Length>> const c = a * b;
  // The numbers below were obtained using Mathematica.
  EXPECT_THAT(
      c.value,
      AlmostEquals(5910974510923886.0 * std::pow(0.5, 48) * Metre * Metre,
                   0));
  EXPECT_THAT(
      c.error,
      AlmostEquals(7881299347837104.0 * std::pow(0.5, 102) * Metre * Metre,
                   0));
}

TEST_F(DoublePrecisionTest, LongQuotient) {
  DoublePrecision<Length> a(3 * Metre);
  a.Increment(474 * ε * Metre);
  DoublePrecision<Length> b(7 * Metre);
  b.Increment(-517 * ε * Metre);
  DoublePrecision<double> const c = a / b;
  // The numbers below were obtained using Mathematica.
  EXPECT_THAT(c.value,
              AlmostEquals(7720456504064105.0 * std::pow(0.5, 54), 0));
  EXPECT_THAT(c.error,
              AlmostEquals(-7352815717686216.0 * std::pow(0.5, 110), 0));
}

TEST_F(DoublePrecisionTest, Mod2π) {
  // Slightly above 2000 π.
  DoublePrecision<Angle> a((2e3 * 103'993 + 2) / 33'102 * Radian);
  auto const c = Mod2π(a);
  // HexLiteral[
  //     CorrectlyRound[Mod[CorrectlyRound[(2000 103993 + 2)/33102], 2 π]]]
  EXPECT_THAT(c.value + c.error,
              AlmostEquals(0x1.F12375A5877D6p-15 * Radian, 0));
  // HexLiteral[CorrectlyRound[(2000 103993 + 2)/33102] -
  //     CorrectlyRound[2000 CorrectlyRound[π]]]
  EXPECT_THAT(a.value - 2000 * π * Radian,
              AlmostEquals(0x1.F123760000000p-15 * Radian, 0));
}

}  // namespace numerics
}  // namespace principia

#undef PRINCIPIA_USE_IACA
