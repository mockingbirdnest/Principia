#include "numerics/wide.hpp"

#include "astronomy/time_scales.hpp"
#include "base/algebra.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "numerics/fma.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using ::testing::Eq;
using ::testing::Ne;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_algebra;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_wide;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;

class WideTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
};

TEST_F(WideTest, Concepts) {
  static_assert(affine_module<Wide<Instant>, Wide<double>>);
  static_assert(ring<Wide<double>>);
  static_assert(module<Displacement<World>, Wide<double>>);
  // The field axioms require 1 / x, Identity() / x is not enough.
  static_assert(!field<Wide<double>>);
}

TEST_F(WideTest, Double) {
  Displacement<World> const r{{2 * Metre, 3 * Metre, 5 * Metre}};
  Wide<double> const ℼ = Wide(π);
  Wide<double> const 𝕖 = Wide(e);
  EXPECT_THAT(static_cast<double>(+𝕖), Eq(e));
  EXPECT_THAT(static_cast<double>(-𝕖), Eq(-e));
  EXPECT_THAT(static_cast<double>(𝕖 + ℼ), Eq(e + π));
  EXPECT_THAT(static_cast<double>(𝕖 - ℼ), Eq(e - π));
  EXPECT_THAT(static_cast<double>(𝕖 * ℼ), Eq(e * π));
  EXPECT_THAT(static_cast<double>(𝕖 / ℼ), Eq(e / π));
  EXPECT_THAT(ℼ * r, Eq(π * r));
  EXPECT_THAT(r * ℼ, Eq(r * π));
  EXPECT_THAT(r / ℼ, Eq(r / π));

  Wide<double> 𝕩 = ℼ;
  Wide<double> 𝕪 = 𝕖;
  𝕩 += 𝕪;
  EXPECT_THAT(static_cast<double>(𝕩), AlmostEquals(e + π, 0));
  𝕩 -= 𝕪;
  EXPECT_THAT(static_cast<double>(𝕩), AlmostEquals(π, 1));
  𝕩 *= 2;
  EXPECT_THAT(static_cast<double>(𝕩), AlmostEquals(2 * π, 1));
  𝕩 *= 𝕪;
  EXPECT_THAT(static_cast<double>(𝕩), AlmostEquals(2 * π * e, 1));
  𝕩 /= 𝕪;
  EXPECT_THAT(static_cast<double>(𝕩), AlmostEquals(2 * π, 2));

  // `Wide` works with `RelativeDegreesOfFreedom` without needing to declare any
  // new operators, because `Displacement` and `Velocity` are modules over
  // `Wide<double>`.
  // Note that the `operator/` requires only a module, see the comments in
  // `direct_sum.hpp`.
  Velocity<World> const v{
      {11 * Metre / Second, 23 * Metre / Second, 29 * Metre / Second}};
  RelativeDegreesOfFreedom<World> rdof = {r, v};
  EXPECT_THAT(ℼ * rdof, Eq(π * rdof));
  EXPECT_THAT(rdof / ℼ, Eq(rdof / π));

  auto w = v;
  w *= 𝕖;
  EXPECT_THAT(w, Eq(e * v));
  w /= 𝕖;
  EXPECT_THAT(w, Eq(v));
}

TEST_F(WideTest, Quantity) {
  Displacement<World> const r{{2 * Metre, 3 * Metre, 5 * Metre}};
  Wide<Mass> const ℼ = Wide(π * Kilogram);
  Wide<Mass> const 𝕖 = Wide(e * Kilogram);
  EXPECT_THAT(static_cast<Mass>(+𝕖), Eq(e * Kilogram));
  EXPECT_THAT(static_cast<Mass>(-𝕖), Eq(-e * Kilogram));
  EXPECT_THAT(static_cast<Mass>(𝕖 + ℼ), Eq((e + π) * Kilogram));
  EXPECT_THAT(static_cast<Mass>(𝕖 - ℼ), Eq((e - π) * Kilogram));
  EXPECT_THAT(static_cast<Square<Mass>>(𝕖 * ℼ), Eq(e * π * Pow<2>(Kilogram)));
  EXPECT_THAT(static_cast<double>(𝕖 / ℼ), Eq(e / π));
  EXPECT_THAT(ℼ * r, Eq(π * Kilogram * r));
  EXPECT_THAT(r * ℼ, Eq(r * π * Kilogram));
  EXPECT_THAT(r / ℼ, Eq(r / (π * Kilogram)));

  Wide<Mass> 𝕩 = ℼ;
  Wide<Mass> 𝕪 = 𝕖;
  𝕩 += 𝕪;
  EXPECT_THAT(static_cast<Mass>(𝕩), AlmostEquals((e + π) * Kilogram, 0));
  𝕩 -= 𝕪;
  EXPECT_THAT(static_cast<Mass>(𝕩), AlmostEquals(π * Kilogram, 1));
  𝕩 *= 2;
  EXPECT_THAT(static_cast<Mass>(𝕩), AlmostEquals(2 * π * Kilogram, 1));

  Velocity<World> const v{
      {11 * Metre / Second, 23 * Metre / Second, 29 * Metre / Second}};
  RelativeDegreesOfFreedom<World> rdof = {r, v};
  EXPECT_THAT(ℼ * rdof, Eq(π * Kilogram * rdof));
  EXPECT_THAT(rdof / ℼ, Eq(rdof / (π * Kilogram)));
}

TEST_F(WideTest, Instant) {
  Instant const t = "2026-10-17T14:47:11"_TT;
  Instant const t₀ = "2026-10-17T00:00:00"_TT;
  Time const Δt = t - t₀;
  Wide<Instant> const 𝕥 = Wide(t);
  Wide<Instant> const 𝕥₀ = Wide(t₀);
  Wide<Time> const Δ𝕥 = Wide(Δt);
  EXPECT_THAT(static_cast<Time>(𝕥 - 𝕥₀), Eq(Δt));
  EXPECT_THAT(static_cast<Instant>(𝕥₀ + Δ𝕥), Eq(t));
  EXPECT_THAT(static_cast<Instant>(𝕥 - Δ𝕥), Eq(t₀));
}

TEST_F(WideTest, FMA) {
  Wide<double> const ℼ = Wide(π);
  Wide<double> const 𝕖 = Wide(e);
  constexpr double x = e + 4 * π;
  constexpr double y = e - π;
  EXPECT_THAT(e * x + π, AlmostEquals(FusedMultiplyAdd(e, x, π), 1));
  EXPECT_THAT(e * y - π, AlmostEquals(FusedMultiplySubtract(e, y, π), 1));
  EXPECT_THAT(-e * x - π,
              AlmostEquals(FusedNegatedMultiplySubtract(e, x, π), 1));
  EXPECT_THAT(-e * y + π, AlmostEquals(FusedNegatedMultiplyAdd(e, y, π), 1));
  Wide<double> const 𝕩 = Wide(x);
  Wide<double> const 𝕪 = Wide(y);
  EXPECT_THAT(static_cast<double>(FusedMultiplyAdd(𝕖, 𝕩, ℼ)),
              Eq(FusedMultiplyAdd(e, x, π)));
  EXPECT_THAT(static_cast<double>(FusedMultiplySubtract(𝕖, 𝕪, ℼ)),
              Eq(FusedMultiplySubtract(e, y, π)));
  EXPECT_THAT(static_cast<double>(FusedNegatedMultiplySubtract(𝕖, 𝕩, ℼ)),
              Eq(FusedNegatedMultiplySubtract(e, x, π)));
  EXPECT_THAT(static_cast<double>(FusedNegatedMultiplyAdd(𝕖, 𝕪, ℼ)),
              Eq(FusedNegatedMultiplyAdd(e, y, π)));

  Displacement<World> const u{{x * Metre, -y * Metre, 3 * x * Metre}};
  Displacement<World> const v{{π * Metre, -π * Metre, 3 * y * Metre}};
  Position<World> const q = World::origin + v;
  EXPECT_THAT(e * u + v, AlmostEquals(FusedMultiplyAdd(e, u, v), 1));
  EXPECT_THAT(e * u - v, AlmostEquals(FusedMultiplySubtract(e, u, v), 1));
  EXPECT_THAT(-e * u - v,
              AlmostEquals(FusedNegatedMultiplySubtract(e, u, v), 1));
  EXPECT_THAT(-e * u + v, AlmostEquals(FusedNegatedMultiplyAdd(e, u, v), 1));

  EXPECT_THAT(FusedMultiplyAdd(𝕖, u, v),
              Eq(FusedMultiplyAdd(e, u, v)));
  EXPECT_THAT(FusedMultiplySubtract(𝕖, u, v),
              Eq(FusedMultiplySubtract(e, u, v)));
  EXPECT_THAT(FusedNegatedMultiplySubtract(𝕖, u, v),
              Eq(FusedNegatedMultiplySubtract(e, u, v)));
  EXPECT_THAT(FusedNegatedMultiplyAdd(𝕖, u, v),
              Eq(FusedNegatedMultiplyAdd(e, u, v)));

  EXPECT_THAT(FusedMultiplyAdd(u, 𝕖, v),
              Eq(FusedMultiplyAdd(u, e, v)));
  EXPECT_THAT(FusedMultiplySubtract(u, 𝕖, v),
              Eq(FusedMultiplySubtract(u, e, v)));
  EXPECT_THAT(FusedNegatedMultiplySubtract(u, 𝕖, v),
              Eq(FusedNegatedMultiplySubtract(u, e, v)));
  EXPECT_THAT(FusedNegatedMultiplyAdd(u, 𝕖, v),
              Eq(FusedNegatedMultiplyAdd(u, e, v)));

  EXPECT_THAT(FusedMultiplyAdd(𝕖, u, q),
              Eq(FusedMultiplyAdd(e, u, q)));
  EXPECT_THAT(FusedNegatedMultiplyAdd(𝕖, u, q),
              Eq(FusedNegatedMultiplyAdd(e, u, q)));

  EXPECT_THAT(FusedMultiplyAdd(u, 𝕖, q),
              Eq(FusedMultiplyAdd(u, e, q)));
  EXPECT_THAT(FusedNegatedMultiplyAdd(u, 𝕖, q),
              Eq(FusedNegatedMultiplyAdd(u, e, q)));
}

}  // namespace numerics
}  // namespace principia
