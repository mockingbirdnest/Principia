#include "numerics/wide.hpp"

#include "astronomy/time_scales.hpp"
#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace numerics {

using ::testing::Eq;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_algebra;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_wide;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_matchers;

class WideTest : public ::testing::Test {
protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
};

TEST_F(WideTest, Concepts) {
  static_assert(ring<Wide<double>>);
  static_assert(module<Displacement<World>, Wide<double>>);
  // The field axioms require 1 / x, Identity() / x is not enough.
  static_assert(!field<Wide<double>>);
}

TEST_F(WideTest, Double) {
  Displacement<World> const r{{2 * Metre, 3 * Metre, 5 * Metre}};
  Wide<double> ℼ = Wide(π);
  Wide<double> 𝕖 = Wide(e);
  EXPECT_THAT((+𝕖).narrow(), Eq(e));
  EXPECT_THAT((-𝕖).narrow(), Eq(-e));
  EXPECT_THAT((𝕖 + ℼ).narrow(), Eq(e+π));
  EXPECT_THAT((𝕖 - ℼ).narrow(), Eq(e-π));
  EXPECT_THAT((𝕖 * ℼ).narrow(), Eq(e*π));
  EXPECT_THAT((𝕖 / ℼ).narrow(), Eq(e/π));
  EXPECT_THAT(ℼ * r, Eq(π * r));
  EXPECT_THAT(r * ℼ, Eq(r * π));
  EXPECT_THAT(r / ℼ, Eq(r / π));

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
}

TEST_F(WideTest, Quantity) {
  Displacement<World> const r{{2 * Metre, 3 * Metre, 5 * Metre}};
  Wide<Mass> ℼ = Wide(π * Kilogram);
  Wide<Mass> 𝕖 = Wide(e * Kilogram);
  EXPECT_THAT((+𝕖).narrow(), Eq(e * Kilogram));
  EXPECT_THAT((-𝕖).narrow(), Eq(-e * Kilogram));
  EXPECT_THAT((𝕖 + ℼ).narrow(), Eq((e + π) * Kilogram));
  EXPECT_THAT((𝕖 - ℼ).narrow(), Eq((e - π) * Kilogram));
  EXPECT_THAT((𝕖 * ℼ).narrow(), Eq(e * π * Kilogram));
  EXPECT_THAT((𝕖 / ℼ).narrow(), Eq(e / π * Kilogram));
  EXPECT_THAT(ℼ * r, Eq(π * r));
  EXPECT_THAT(r * ℼ, Eq(r * π));
  EXPECT_THAT(r / ℼ, Eq(r / π));

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
  EXPECT_THAT((𝕥 - 𝕥₀).narrow(), Eq(Δt));
  EXPECT_THAT((𝕥₀ + Δ𝕥).narrow(), Eq(t));
  EXPECT_THAT((𝕥 - Δ𝕥).narrow(), Eq(t₀));
}

}  // namespace quantities
}  // namespace principia
