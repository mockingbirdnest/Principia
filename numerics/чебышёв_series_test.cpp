#include "numerics/чебышёв_series.hpp"

#include "astronomy/frames.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/matrix_computations.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {

using ::testing::ElementsAre;
using ::testing::Lt;
using namespace principia::astronomy::_frames;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::numerics::_чебышёв_series;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_numerics_matchers;

class ЧебышёвSeriesTest : public ::testing::Test {
 protected:
  ЧебышёвSeriesTest()
      : t_min_(t0_ - 1 * Second),
        t_max_(t0_ + 3 * Second) {}

  Instant const t0_;
  Instant t_min_;
  Instant t_max_;
};

using ЧебышёвSeriesDeathTest = ЧебышёвSeriesTest;

TEST_F(ЧебышёвSeriesDeathTest, ConstructionErrors) {
  using Series = ЧебышёвSeries<double, Instant>;
  EXPECT_DEATH({ Series p({}, t_min_, t_max_); }, "at least 0");
  EXPECT_DEATH({ Series p({1}, t_max_, t_min_); }, "not be empty");
}

TEST_F(ЧебышёвSeriesDeathTest, EvaluationErrors) {
#ifdef _DEBUG
  using Series = ЧебышёвSeries<double, Instant>;
  EXPECT_DEATH({
    Series p({1}, t_min_, t_max_);
    p.Evaluate(t_min_ - 10 * Second);
  }, ">= -1.1");
  EXPECT_DEATH({
    Series p({1}, t_min_, t_max_);
    p.Evaluate(t_max_ + 10 * Second);
  }, "<= 1.1");
#endif
}

TEST_F(ЧебышёвSeriesTest, T0) {
  ЧебышёвSeries<double, Instant> t0({1}, t_min_, t_max_);
  EXPECT_EQ(1, t0.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(1, t0.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, T1) {
  ЧебышёвSeries<double, Instant> t1({0, 1}, t_min_, t_max_);
  EXPECT_EQ(0, t1.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(1, t1.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, T2) {
  ЧебышёвSeries<double, Instant> t2({0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(1, t2.Evaluate(t0_ + -1 * Second));
  EXPECT_EQ(-1, t2.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(1, t2.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, T3) {
  ЧебышёвSeries<double, Instant> t3({0, 0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(-1, t3.Evaluate(t0_ + -1 * Second));
  EXPECT_EQ(0, t3.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(-1, t3.Evaluate(t0_ + 2 * Second));
  EXPECT_EQ(1, t3.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, X5) {
  ЧебышёвSeries<double, Instant> x5(
      {0.0, 10.0 / 16.0, 0, 5.0 / 16.0, 0, 1.0 / 16.0},
      t_min_, t_max_);
  EXPECT_EQ(-1, x5.Evaluate(t0_ + -1 * Second));
  EXPECT_EQ(0, x5.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(1.0 / 1024.0, x5.Evaluate(t0_ + 1.5 * Second));
  EXPECT_EQ(1.0 / 32.0, x5.Evaluate(t0_ + 2 * Second));
  EXPECT_EQ(1, x5.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, X6) {
  ЧебышёвSeries<double, Instant> x6(
      {10.0 / 32.0, 0, 15.0 / 32.0, 0, 6.0 / 32.0, 0, 1.0 / 32.0},
      t_min_, t_max_);
  EXPECT_EQ(1, x6.Evaluate(t0_ + -1 * Second));
  EXPECT_EQ(0, x6.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(1.0 / 4096.0, x6.Evaluate(t0_ + 1.5 * Second));
  EXPECT_EQ(1.0 / 64.0, x6.Evaluate(t0_ + 2 * Second));
  EXPECT_EQ(1, x6.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, T2Dimension) {
  ЧебышёвSeries<Length, Instant> t2(
      {0 * Metre, 0 * Metre, 1 * Metre},
      t_min_, t_max_);
  EXPECT_EQ(1 * Metre, t2.Evaluate(t0_ + -1 * Second));
  EXPECT_EQ(-1 * Metre, t2.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(1 * Metre, t2.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesTest, T2Double) {
  ЧебышёвSeries<Length, double> t2(
      {0 * Metre, 0 * Metre, 1 * Metre},
      -1, 2);
  EXPECT_EQ(1 * Metre, t2.Evaluate(-1));
  EXPECT_EQ(-7.0 / 9.0 * Metre, t2.Evaluate(1));
  EXPECT_EQ(1 * Metre, t2.Evaluate(2));
}

TEST_F(ЧебышёвSeriesTest, FrobeniusCompanionMatrix) {
  ЧебышёвSeries<double, Instant> series({-2, 3, 5, 6}, t_min_, t_max_);
  auto const matrix = series.FrobeniusCompanionMatrix();
  EXPECT_THAT(matrix,
              AlmostEquals(UnboundedMatrix<double>(
                  {0.0,             1.0,         0.0,
                   1.0 / 2.0,       0.0,   1.0 / 2.0,
                   1.0 / 6.0, 1.0 / 4.0, -5.0 / 12.0}),
                  0));
  auto const matrix_schur_decomposition = RealSchurDecomposition(matrix, 1e-16);
  EXPECT_THAT(matrix_schur_decomposition.real_eigenvalues,
              ElementsAre(AlmostEquals((1.0 - Sqrt(337.0)) / 24.0, 4),
                          AlmostEquals(-0.5, 1),
                          AlmostEquals((1.0 + Sqrt(337.0)) / 24.0, 2)));
}

TEST_F(ЧебышёвSeriesTest, MayHaveRealRoots) {
  // B₀ path.
  ЧебышёвSeries<double, Instant> series1({16, 5, 3, 7}, t_min_, t_max_);
  EXPECT_FALSE(series1.MayHaveRealRoots());
  // We don't know, but it actually doesn't have zeroes.
  ЧебышёвSeries<double, Instant> series3({13, 5, 3, 7}, t_min_, t_max_);
  EXPECT_TRUE(series3.MayHaveRealRoots());
  // We don't know, but it actually has zeroes.
  ЧебышёвSeries<double, Instant> series4({4, 5, 3, 7}, t_min_, t_max_);
  EXPECT_TRUE(series4.MayHaveRealRoots());
}

TEST_F(ЧебышёвSeriesTest, RealRoots) {
  ЧебышёвSeries<double, Instant> series({-2, 3, 5, 6}, t_min_, t_max_);
  Instant const r1 = t0_ + (13.0 - Sqrt(337.0)) / 12.0 * Second;
  Instant const r2 = t0_;
  Instant const r3 = t0_ + (13.0 + Sqrt(337.0)) / 12.0 * Second;
  EXPECT_THAT(series.RealRoots(1e-16),
              ElementsAre(AlmostEquals(r1, 15),
                          AbsoluteErrorFrom(r2, Lt(2.3e-16 * Second)),
                          AlmostEquals(r3, 1)));
  EXPECT_THAT(Abs(series.Evaluate(r1)), Lt(8.9e-16));
  EXPECT_THAT(Abs(series.Evaluate(r2)), AlmostEquals(0, 0));
  EXPECT_THAT(Abs(series.Evaluate(r3)), Lt(1.8e-15));
}

TEST_F(ЧебышёвSeriesTest, X6Vector) {
  using V = Vector<Length, ICRS>;
  // {T3, X5, X6}
  V const c0 = V({0.0 * Metre, 0.0 * Metre, 10.0 / 32.0 * Metre});
  V const c1 = V({0.0 * Metre, 10.0 / 16.0 * Metre, 0.0 * Metre});
  V const c2 = V({0.0 * Metre, 0.0 * Metre, 15.0 / 32.0 * Metre});
  V const c3 = V({1.0 * Metre, 5.0 / 16.0 * Metre, 0.0 * Metre});
  V const c4 = V({0.0 * Metre, 0.0 * Metre, 6.0 / 32.0 * Metre});
  V const c5 = V({0.0 * Metre, 1.0 / 16.0 * Metre, 0 * Metre});
  V const c6 = V({0.0 * Metre, 0.0 * Metre, 1.0 / 32.0 * Metre});
  ЧебышёвSeries<Vector<Length, ICRS>, Instant> x6(
      {c0, c1, c2, c3, c4, c5, c6},
      t_min_, t_max_);
  EXPECT_EQ(V({-1 * Metre, -1 * Metre, 1 * Metre}),
            x6.Evaluate(t0_ + -1 * Second));
  EXPECT_EQ(V({0 * Metre, 0 * Metre, 0 * Metre}),
            x6.Evaluate(t0_ + 1 * Second));
  EXPECT_EQ(V({-1 * Metre, 1.0 / 32.0 * Metre, 1 / 64.0 * Metre}),
            x6.Evaluate(t0_ + 2 * Second));
  EXPECT_EQ(V({1 * Metre, 1 * Metre, 1 * Metre}),
            x6.Evaluate(t0_ + 3 * Second));
}

TEST_F(ЧебышёвSeriesDeathTest, SerializationError) {
  using Series1 = ЧебышёвSeries<Speed, Instant>;
  using Series2 = ЧебышёвSeries<double, Instant>;
  Series1 v({1 * Metre / Second,
             -2 * Metre / Second,
             5 * Metre / Second},
             t_min_, t_max_);
  Series2 d({7, 8, -1}, t_min_, t_max_);

  EXPECT_DEATH({
    serialization::ЧебышёвSeries message;
    v.WriteToMessage(&message);
    Series2::ReadFromMessage(message);
  }, "has_double");
  EXPECT_DEATH({
    serialization::ЧебышёвSeries message;
    d.WriteToMessage(&message);
    Series1::ReadFromMessage(message);
  }, "has_quantity");
}

TEST_F(ЧебышёвSeriesTest, SerializationSuccess) {
  {
    serialization::ЧебышёвSeries message;
    ЧебышёвSeries<Speed, Instant> const v1({1 * Metre / Second,
                                            -2 * Metre / Second,
                                            5 * Metre / Second},
                                           t_min_, t_max_);
    v1.WriteToMessage(&message);
    EXPECT_EQ(3, message.coefficient_size());
    EXPECT_FALSE(message.coefficient(0).has_double_());
    EXPECT_TRUE(message.coefficient(0).has_quantity());
    EXPECT_EQ(0x7C01, message.coefficient(0).quantity().dimensions());
    EXPECT_EQ(1.0, message.coefficient(0).quantity().magnitude());
    EXPECT_TRUE(message.has_lower_bound());
    EXPECT_TRUE(message.lower_bound().has_scalar());
    EXPECT_TRUE(message.lower_bound().scalar().has_dimensions());
    EXPECT_TRUE(message.lower_bound().scalar().has_magnitude());
    EXPECT_EQ(-1.0, message.lower_bound().scalar().magnitude());
    EXPECT_TRUE(message.has_upper_bound());
    EXPECT_TRUE(message.upper_bound().has_scalar());
    EXPECT_TRUE(message.upper_bound().scalar().has_dimensions());
    EXPECT_TRUE(message.upper_bound().scalar().has_magnitude());
    EXPECT_EQ(3.0, message.upper_bound().scalar().magnitude());
    auto const v2 = ЧебышёвSeries<Speed, Instant>::ReadFromMessage(message);
    EXPECT_EQ(v1, v2);
  }
  {
    serialization::ЧебышёвSeries message;
    ЧебышёвSeries<double, Instant> const d1({-1, 2, 5}, t_min_, t_max_);
    d1.WriteToMessage(&message);
    EXPECT_EQ(3, message.coefficient_size());
    EXPECT_TRUE(message.coefficient(0).has_double_());
    EXPECT_FALSE(message.coefficient(0).has_quantity());
    EXPECT_EQ(-1.0, message.coefficient(0).double_());
    EXPECT_TRUE(message.has_lower_bound());
    EXPECT_TRUE(message.lower_bound().has_scalar());
    EXPECT_TRUE(message.lower_bound().scalar().has_dimensions());
    EXPECT_TRUE(message.lower_bound().scalar().has_magnitude());
    EXPECT_EQ(-1.0, message.lower_bound().scalar().magnitude());
    EXPECT_TRUE(message.has_upper_bound());
    EXPECT_TRUE(message.upper_bound().has_scalar());
    EXPECT_TRUE(message.upper_bound().scalar().has_dimensions());
    EXPECT_TRUE(message.upper_bound().scalar().has_magnitude());
    EXPECT_EQ(3.0, message.upper_bound().scalar().magnitude());
    auto const d2 = ЧебышёвSeries<double, Instant>::ReadFromMessage(message);
    EXPECT_EQ(d1, d2);
  }
}

}  // namespace numerics
}  // namespace principia
