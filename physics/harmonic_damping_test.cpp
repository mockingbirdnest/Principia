#include "physics/harmonic_damping.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Vector;
using quantities::Inverse;
using quantities::Length;
using quantities::Pow;
using quantities::Square;
using quantities::si::Metre;
using ::testing::Eq;

class HarmonicDampingTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
};

TEST_F(HarmonicDampingTest, AccelerationQuantities) {
  HarmonicDamping σ(1 * Metre);
  EXPECT_THAT(σ.inner_threshold(), Eq(1 * Metre));
  EXPECT_THAT(σ.outer_threshold(), Eq(3 * Metre));
  Vector<double, World> x({1, 0, 0});
  Inverse<Square<Length>> const ℜ_over_r = 5 / Pow<2>(Metre);
  Inverse<Square<Length>> const ℜʹ = 17 / Pow<2>(Metre);
  Inverse<Square<Length>> σℜ_over_r;
  Vector<Inverse<Square<Length>>, World> grad_σℜ;

  {
    Length const r = 3 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    EXPECT_THAT(σℜ_over_r, Eq(0 / Pow<2>(Metre)));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(0 / Pow<2>(Metre)));
  }
  {
    Length const r = 2 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    auto const ℜ = ℜ_over_r * r;
    auto const σ = 0.5;
    auto const σʹ = -3 / (4 * Metre);
    EXPECT_THAT(σℜ_over_r, Eq(ℜ_over_r / 2));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(σʹ * ℜ + ℜʹ * σ));
  }
  {
    Length const r = 1 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    EXPECT_THAT(σℜ_over_r, Eq(ℜ_over_r));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(ℜʹ));
  }
  {
    Length const r = 0.5 * Metre;
    σ.ComputeDampedRadialQuantities(
        r, r * r, x, ℜ_over_r, ℜʹ, σℜ_over_r, grad_σℜ);
    EXPECT_THAT(σℜ_over_r, Eq(ℜ_over_r));
    EXPECT_THAT(grad_σℜ.coordinates().x, Eq(ℜʹ));
  }
}

}  // namespace physics
}  // namespace principia
