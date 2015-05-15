
#include "numerics/чебышёв_series.hpp"

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using si::Second;

namespace numerics {

class ЧебышёвSeriesTest : public ::testing::Test {
 protected:
  ЧебышёвSeriesTest()
      : t_min_(-1 * Second),
        t_max_(3 * Second) {}

  Instant t_min_;
  Instant t_max_;
};

TEST_F(ЧебышёвSeriesTest, T0) {
  ЧебышёвSeries<double> t0({1}, t_min_, t_max_);
  EXPECT_EQ(1, t0.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1, t0.Evaluate(Instant(3 * Second)));
}

TEST_F(ЧебышёвSeriesTest, T1) {
  ЧебышёвSeries<double> t1({0, 1}, t_min_, t_max_);
  EXPECT_EQ(0, t1.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1, t1.Evaluate(Instant(3 * Second)));
}

TEST_F(ЧебышёвSeriesTest, T2) {
  ЧебышёвSeries<double> t2({0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(1, t2.Evaluate(Instant(-1 * Second)));
  EXPECT_EQ(-1, t2.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1, t2.Evaluate(Instant(3 * Second)));
}

TEST_F(ЧебышёвSeriesTest, T3) {
  ЧебышёвSeries<double> t3({0, 0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(-1, t3.Evaluate(Instant(-1 * Second)));
  EXPECT_EQ(0, t3.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(-1, t3.Evaluate(Instant(2 * Second)));
  EXPECT_EQ(1, t3.Evaluate(Instant(3 * Second)));
}

TEST_F(ЧебышёвSeriesTest, X5) {
  ЧебышёвSeries<double> x5({0.0, 10.0 / 16.0, 0, 5.0 / 16.0, 0, 1.0 / 16.0},
                           t_min_, t_max_);
  EXPECT_EQ(-1, x5.Evaluate(Instant(-1 * Second)));
  EXPECT_EQ(0, x5.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1.0 / 1024.0, x5.Evaluate(Instant(1.5 * Second)));
  EXPECT_EQ(1.0 / 32.0, x5.Evaluate(Instant(2 * Second)));
  EXPECT_EQ(1, x5.Evaluate(Instant(3 * Second)));
}

TEST_F(ЧебышёвSeriesTest, X6) {
  ЧебышёвSeries<double> x6(
      {10.0 / 32.0, 0, 15.0 / 32.0, 0, 6.0 / 32.0, 0, 1.0 / 32.0},
      t_min_, t_max_);
  EXPECT_EQ(1, x6.Evaluate(Instant(-1 * Second)));
  EXPECT_EQ(0, x6.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1.0 / 4096.0, x6.Evaluate(Instant(1.5 * Second)));
  EXPECT_EQ(1.0 / 64.0, x6.Evaluate(Instant(2 * Second)));
  EXPECT_EQ(1, x6.Evaluate(Instant(3 * Second)));
}

}  // namespace numerics
}  // namespace principia
