
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
using quantities::Temperature;
using quantities::si::Kelvin;
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

  using Series2 = PoissonSeries<Temperature, 2, 2, HornerEvaluator>;
  using Series3 = PoissonSeries<Displacement<World>, 3, 3, HornerEvaluator>;

  Instant const t0_;
};

TEST_F(PoissonSeriesBasisTest, AperiodicScalar) {
  auto const aperiodic = PoissonSeriesBasisGenerator<
      Series2,
      /*degree=*/2>::Basis(t0_);
  EXPECT_EQ(3, aperiodic.size());

  Instant const t1 = t0_ + 2 * Second;

  EXPECT_EQ(1 * Kelvin, aperiodic[0](t1));
  EXPECT_EQ(2 * Kelvin, aperiodic[1](t1));
  EXPECT_EQ(4 * Kelvin, aperiodic[2](t1));
}

TEST_F(PoissonSeriesBasisTest, AperiodicVector) {
  auto const aperiodic = PoissonSeriesBasisGenerator<
      Series3,
      /*degree=*/3>::Basis(t0_);
  auto const aperiodic_subspaces =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/3>::Subspaces(t0_);
  EXPECT_EQ(12, aperiodic.size());

  Instant const t1 = t0_ + 2 * Second;

  // Degree 0.
  EXPECT_EQ(Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[0](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre}),
            aperiodic[1](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre}),
            aperiodic[2](t1));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                 aperiodic_subspaces[0]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[1]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[2]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                 aperiodic_subspaces[1]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[2]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                 aperiodic_subspaces[2]));

  // Degree 1.
  EXPECT_EQ(Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[3](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
            aperiodic[4](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
            aperiodic[5](t1));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[5]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[5]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[5]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                 aperiodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[5]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                 aperiodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                aperiodic_subspaces[5]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                 aperiodic_subspaces[5]));

  // Degree 2.
  EXPECT_EQ(Displacement<World>({4 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[6](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 4 * Metre, 0 * Metre}),
            aperiodic[7](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 4 * Metre}),
            aperiodic[8](t1));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                 aperiodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[6]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                 aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[7]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                 aperiodic_subspaces[8]));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                aperiodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                aperiodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                aperiodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                aperiodic_subspaces[8]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[6],
                                                 aperiodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[6],
                                                aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[6],
                                                aperiodic_subspaces[8]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[7],
                                                 aperiodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[7],
                                                aperiodic_subspaces[8]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[8],
                                                 aperiodic_subspaces[8]));

  // Degree 3.
  EXPECT_EQ(Displacement<World>({8 * Metre, 0 * Metre, 0 * Metre}),
            aperiodic[9](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 8 * Metre, 0 * Metre}),
            aperiodic[10](t1));
  EXPECT_EQ(Displacement<World>({0 * Metre, 0 * Metre, 8 * Metre}),
            aperiodic[11](t1));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[0],
                                                aperiodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[1],
                                                aperiodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[2],
                                                aperiodic_subspaces[11]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                 aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[3],
                                                aperiodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                aperiodic_subspaces[9]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                 aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[4],
                                                aperiodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                aperiodic_subspaces[10]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[5],
                                                 aperiodic_subspaces[11]));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[6],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[6],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[6],
                                                aperiodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[7],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[7],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[7],
                                                aperiodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[8],
                                                aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[8],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[8],
                                                aperiodic_subspaces[11]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[9],
                                                 aperiodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[9],
                                                aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[9],
                                                aperiodic_subspaces[11]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[10],
                                                 aperiodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[10],
                                                aperiodic_subspaces[11]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(aperiodic_subspaces[11],
                                                 aperiodic_subspaces[11]));
}

TEST_F(PoissonSeriesBasisTest, PeriodicScalar) {
  AngularFrequency const ω = π / 6 * Radian / Second;
  auto const periodic = PoissonSeriesBasisGenerator<
      Series2,
      /*degree=*/2>::Basis(ω, t0_);
  EXPECT_EQ(6, periodic.size());

  Instant const t1 = t0_ + 2 * Second;

  EXPECT_THAT(periodic[0](t1), AlmostEquals(0.5 * Kelvin, 1));
  EXPECT_THAT(periodic[1](t1), AlmostEquals(Sqrt(3) / 2 * Kelvin, 0));

  EXPECT_THAT(periodic[2](t1), AlmostEquals(1 * Kelvin, 1));
  EXPECT_THAT(periodic[3](t1), AlmostEquals(Sqrt(3) * Kelvin, 0));

  EXPECT_THAT(periodic[4](t1), AlmostEquals(2 * Kelvin, 1));
  EXPECT_THAT(periodic[5](t1), AlmostEquals(2 * Sqrt(3) * Kelvin, 0));
}

TEST_F(PoissonSeriesBasisTest, PeriodicVector) {
  AngularFrequency const ω = π / 6 * Radian / Second;
  auto const periodic = PoissonSeriesBasisGenerator<
      Series3,
      /*degree=*/3>::Basis(ω, t0_);
  auto const periodic_subspaces =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/3>::Subspaces(ω, t0_);
  EXPECT_EQ(24, periodic.size());

  Instant const t1 = t0_ + 2 * Second;

  // Degree 0, Cos.
  EXPECT_THAT(
      periodic[0](t1),
      AlmostEquals(
          Displacement<World>({0.5 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[1](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0.5 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[2](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 0.5 * Metre}), 1));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                 periodic_subspaces[0]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[1]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[2]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                 periodic_subspaces[1]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[2]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                 periodic_subspaces[2]));

  // Degree 0, Sin.
  EXPECT_THAT(
      periodic[3](t1),
      AlmostEquals(
          Displacement<World>({Sqrt(3) / 2 * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[4](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, Sqrt(3) / 2 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[5](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, Sqrt(3) / 2 * Metre}), 0));

  LOG(ERROR)<<periodic_subspaces[1];
  LOG(ERROR)<<periodic_subspaces[4];
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[5]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[5]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[5]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                 periodic_subspaces[3]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[5]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                 periodic_subspaces[4]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                periodic_subspaces[5]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                 periodic_subspaces[5]));

  // Degree 1, Cos.
  EXPECT_THAT(
      periodic[6](t1),
      AlmostEquals(
          Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[7](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[8](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre}), 1));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[8]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                 periodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                periodic_subspaces[6]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                 periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                periodic_subspaces[8]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                periodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                periodic_subspaces[7]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                 periodic_subspaces[8]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[6],
                                                 periodic_subspaces[6]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[6],
                                                periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[6],
                                                periodic_subspaces[8]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[7],
                                                 periodic_subspaces[7]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[7],
                                                periodic_subspaces[8]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[8],
                                                 periodic_subspaces[8]));

  // Degree 1, Sin.
  EXPECT_THAT(
      periodic[9](t1),
      AlmostEquals(
          Displacement<World>({Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[10](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[11](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, Sqrt(3) * Metre}), 0));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                 periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[0],
                                                periodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[9]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                 periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[1],
                                                periodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                periodic_subspaces[10]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[2],
                                                 periodic_subspaces[11]));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[3],
                                                periodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[4],
                                                periodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[5],
                                                periodic_subspaces[11]));

  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[6],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[6],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[6],
                                                periodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[7],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[7],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[7],
                                                periodic_subspaces[11]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[8],
                                                periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[8],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[8],
                                                periodic_subspaces[11]));

  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[9],
                                                 periodic_subspaces[9]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[9],
                                                periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[9],
                                                periodic_subspaces[11]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[10],
                                                 periodic_subspaces[10]));
  EXPECT_TRUE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[10],
                                                periodic_subspaces[11]));
  EXPECT_FALSE(PoissonSeriesSubspace::orthogonal(periodic_subspaces[11],
                                                 periodic_subspaces[11]));

  // Degree 2, Cos.
  EXPECT_THAT(
      periodic[12](t1),
      AlmostEquals(
          Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[13](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[14](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}), 1));

  // Degree 2, Sin.
  EXPECT_THAT(
      periodic[15](t1),
      AlmostEquals(
          Displacement<World>({2 * Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[16](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 2 * Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[17](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Sqrt(3) * Metre}), 0));

  // Degree 3, Cos.
  EXPECT_THAT(
      periodic[18](t1),
      AlmostEquals(
          Displacement<World>({4 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[19](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 4 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[20](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 4 * Metre}), 1));

  // Degree 3, Sin.
  EXPECT_THAT(
      periodic[21](t1),
      AlmostEquals(
          Displacement<World>({4 * Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[22](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 4 * Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[23](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 4 * Sqrt(3) * Metre}), 0));
}

TEST_F(PoissonSeriesBasisTest, ReducedDegree) {
  auto const aperiodic = PoissonSeriesBasisGenerator<
      Series3,
      /*degree=*/2>::Basis(t0_);
  EXPECT_EQ(9, aperiodic.size());

  AngularFrequency const ω = π / 6 * Radian / Second;
  auto const periodic = PoissonSeriesBasisGenerator<
      Series3,
      /*degree=*/2>::Basis(ω, t0_);
  EXPECT_EQ(18, periodic.size());

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

  EXPECT_THAT(
      periodic[0](t1),
      AlmostEquals(
          Displacement<World>({0.5 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[1](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0.5 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[2](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 0.5 * Metre}), 1));
  EXPECT_THAT(
      periodic[3](t1),
      AlmostEquals(
          Displacement<World>({Sqrt(3) / 2 * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[4](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, Sqrt(3) / 2 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[5](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, Sqrt(3) / 2 * Metre}), 0));

  EXPECT_THAT(
      periodic[6](t1),
      AlmostEquals(
          Displacement<World>({1 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[7](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 1 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[8](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 1 * Metre}), 1));
  EXPECT_THAT(
      periodic[9](t1),
      AlmostEquals(
          Displacement<World>({Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[10](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[11](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, Sqrt(3) * Metre}), 0));

  EXPECT_THAT(
      periodic[12](t1),
      AlmostEquals(
          Displacement<World>({2 * Metre, 0 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[13](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}), 1));
  EXPECT_THAT(
      periodic[14](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}), 1));
  EXPECT_THAT(
      periodic[15](t1),
      AlmostEquals(
          Displacement<World>({2 * Sqrt(3) * Metre, 0 * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[16](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 2 * Sqrt(3) * Metre, 0 * Metre}), 0));
  EXPECT_THAT(
      periodic[17](t1),
      AlmostEquals(
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Sqrt(3) * Metre}), 0));
}

}  // namespace numerics
}  // namespace principia
