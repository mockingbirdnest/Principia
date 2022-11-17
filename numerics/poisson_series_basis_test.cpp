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

using geometry::Frame;
using geometry::Handedness;
using geometry::Hilbert;
using geometry::Inertial;
using geometry::Instant;
using geometry::Vector;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Sqrt;
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

  using V = Vector<double, World>;
  using Series2 = PoissonSeries<double, 2, 2, HornerEvaluator>;
  using Series3 = PoissonSeries<V, 3, 3, HornerEvaluator>;

  Instant const t0_;
  Instant const t1_ = t0_ - 2 * Second;
  Instant const t2_ = t0_ + 2 * Second;
};

TEST_F(PoissonSeriesBasisTest, AperiodicScalar) {
  auto const aperiodic =
      PoissonSeriesBasisGenerator<Series2,
                                  /*degree=*/2>::Basis(t1_, t2_);
  EXPECT_EQ(3, aperiodic.size());

  EXPECT_EQ(1, aperiodic[0](t2_));
  EXPECT_EQ(1, aperiodic[1](t2_));
  EXPECT_EQ(1, aperiodic[2](t2_));
}

TEST_F(PoissonSeriesBasisTest, AperiodicVector) {
  auto const aperiodic =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/3>::Basis(t1_, t2_);
  auto const aperiodic_subspaces =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/3>::Subspaces();
  EXPECT_EQ(12, aperiodic.size());

  // Degree 0.
  EXPECT_EQ(V({1, 0, 0}), aperiodic[0](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[1](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[2](t2_));

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
  EXPECT_EQ(V({1, 0, 0}), aperiodic[3](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[4](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[5](t2_));

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
  EXPECT_EQ(V({1, 0, 0}), aperiodic[6](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[7](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[8](t2_));

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
  EXPECT_EQ(V({1, 0, 0}), aperiodic[9](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[10](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[11](t2_));

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
  auto const periodic =
      PoissonSeriesBasisGenerator<Series2,
                                  /*degree=*/2>::Basis(ω, t1_, t2_);
  EXPECT_EQ(6, periodic.size());

  EXPECT_THAT(periodic[0](t2_), AlmostEquals(0.5, 1));
  EXPECT_THAT(periodic[1](t2_), AlmostEquals(Sqrt(3) / 2, 0));

  EXPECT_THAT(periodic[2](t2_), AlmostEquals(0.5, 1));
  EXPECT_THAT(periodic[3](t2_), AlmostEquals(Sqrt(3) / 2, 0));

  EXPECT_THAT(periodic[4](t2_), AlmostEquals(0.5, 1));
  EXPECT_THAT(periodic[5](t2_), AlmostEquals(Sqrt(3) / 2, 0));
}

TEST_F(PoissonSeriesBasisTest, PeriodicVector) {
  AngularFrequency const ω = π / 6 * Radian / Second;
  auto const periodic =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/3>::Basis(ω, t1_, t2_);
  auto const periodic_subspaces =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/3>::Subspaces(ω);
  EXPECT_EQ(24, periodic.size());

  // Degree 0, Cos.
  EXPECT_THAT(periodic[0](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[1](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[2](t2_), AlmostEquals(V({0, 0, 0.5}), 1));

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
  EXPECT_THAT(periodic[3](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[4](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[5](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));

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
  EXPECT_THAT(periodic[6](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[7](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[8](t2_), AlmostEquals(V({0, 0, 0.5}), 1));

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
  EXPECT_THAT(periodic[9](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[10](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[11](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));

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
  EXPECT_THAT(periodic[12](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[13](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[14](t2_), AlmostEquals(V({0, 0, 0.5}), 1));

  // Degree 2, Sin.
  EXPECT_THAT(periodic[15](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[16](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[17](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));

  // Degree 3, Cos.
  EXPECT_THAT(periodic[18](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[19](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[20](t2_), AlmostEquals(V({0, 0, 0.5}), 1));

  // Degree 3, Sin.
  EXPECT_THAT(periodic[21](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[22](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[23](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));
}

TEST_F(PoissonSeriesBasisTest, ReducedDegree) {
  auto const aperiodic =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/2>::Basis(t1_, t2_);
  EXPECT_EQ(9, aperiodic.size());

  AngularFrequency const ω = π / 6 * Radian / Second;
  auto const periodic =
      PoissonSeriesBasisGenerator<Series3,
                                  /*degree=*/2>::Basis(ω, t1_, t2_);
  EXPECT_EQ(18, periodic.size());

  EXPECT_EQ(V({1, 0, 0}), aperiodic[0](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[1](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[2](t2_));

  EXPECT_EQ(V({1, 0, 0}), aperiodic[3](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[4](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[5](t2_));

  EXPECT_EQ(V({1, 0, 0}), aperiodic[6](t2_));
  EXPECT_EQ(V({0, 1, 0}), aperiodic[7](t2_));
  EXPECT_EQ(V({0, 0, 1}), aperiodic[8](t2_));

  EXPECT_THAT(periodic[0](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[1](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[2](t2_), AlmostEquals(V({0, 0, 0.5}), 1));
  EXPECT_THAT(periodic[3](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[4](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[5](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));

  EXPECT_THAT(periodic[6](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[7](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[8](t2_), AlmostEquals(V({0, 0, 0.5}), 1));
  EXPECT_THAT(periodic[9](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[10](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[11](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));

  EXPECT_THAT(periodic[12](t2_), AlmostEquals(V({0.5, 0, 0}), 1));
  EXPECT_THAT(periodic[13](t2_), AlmostEquals(V({0, 0.5, 0}), 1));
  EXPECT_THAT(periodic[14](t2_), AlmostEquals(V({0, 0, 0.5}), 1));
  EXPECT_THAT(periodic[15](t2_), AlmostEquals(V({Sqrt(3) / 2, 0, 0}), 0));
  EXPECT_THAT(periodic[16](t2_), AlmostEquals(V({0, Sqrt(3) / 2, 0}), 0));
  EXPECT_THAT(periodic[17](t2_), AlmostEquals(V({0, 0, Sqrt(3) / 2}), 0));
}

}  // namespace numerics
}  // namespace principia
