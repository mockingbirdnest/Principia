
#include "numerics/poisson_series_basis.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Hilbert;
using geometry::Inertial;
using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;

class PoissonSeriesBasisTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Series3 = PoissonSeries<Displacement<World>, 3, HornerEvaluator>;

  Instant const t0_;
};

TEST_F(PoissonSeriesBasisTest, Aperiodic) {
  auto const aperiodic = PoissonSeriesBasisGenerator<
      Series3,
      Hilbert<Displacement<World>>::dimension>::Basis(t0_);
  EXPECT_EQ(12, aperiodic.size());

  Instant const t1 = t0_ + 2 * Second;

  EXPECT_EQ(Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[0](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre}),
            aperiodic[1](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre}),
            aperiodic[2](t1));

  EXPECT_EQ(Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[3](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
            aperiodic[4](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
            aperiodic[5](t1));

  EXPECT_EQ(Displacement<World>({4 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[6](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 4 * Metre, 0 * Metre}),
            aperiodic[7](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 4 * Metre}),
            aperiodic[8](t1));

  EXPECT_EQ(Displacement<World>({8 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[9](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 8 * Metre, 0 * Metre}),
            aperiodic[10](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 8 * Metre}),
            aperiodic[11](t1));
}

TEST_F(PoissonSeriesBasisTest, Periodic) {
  AngularFrequency const ω = π / 6 * Radian / Second;
  auto const periodic = PoissonSeriesBasisGenerator<
      Series3,
      Hilbert<Displacement<World>>::dimension>::Basis(ω, t0_);
  EXPECT_EQ(24, periodic.size());

  Instant const t1 = t0_ + 2 * Second;

  EXPECT_THAT(
      periodic[0](t1),
      AlmostEquals(
          Displacement<World>({Sqrt(3) / 2 * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[1](t1),
      AlmostEquals(
          Displacement<World>({0.5 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[2](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, Sqrt(3) / 2 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[3](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0.5 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[4](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, Sqrt(3) / 2 * Metre}), 0));
  EXPECT_THAT(
      periodic[5](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 0.5 * Metre}), 1));

  EXPECT_THAT(
      periodic[6](t1),
      AlmostEquals(
          Displacement<World>({Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[7](t1),
      AlmostEquals(
          Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[8](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[9](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[10](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, Sqrt(3) * Metre}), 0));
  EXPECT_THAT(
      periodic[11](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre}), 1));

  EXPECT_THAT(
      periodic[12](t1),
      AlmostEquals(
          Displacement<World>({2 * Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[13](t1),
      AlmostEquals(
          Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[14](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 2 * Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[15](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[16](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Sqrt(3) * Metre}), 0));
  EXPECT_THAT(
      periodic[17](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}), 1));

  EXPECT_THAT(
      periodic[18](t1),
      AlmostEquals(
          Displacement<World>({4 * Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[19](t1),
      AlmostEquals(
          Displacement<World>({4 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[20](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 4 * Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[21](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 4 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[22](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 4 * Sqrt(3) * Metre}), 0));
  EXPECT_THAT(
      periodic[23](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 4 * Metre}), 1));
}

}  // namespace numerics
}  // namespace principia
