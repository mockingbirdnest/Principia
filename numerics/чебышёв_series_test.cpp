
#include "numerics/чебышёв_series.hpp"

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using quantities::Speed;
using si::Metre;
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

using ЧебышёвSeriesDeathTest = ЧебышёвSeriesTest;

TEST_F(ЧебышёвSeriesTest, ConstructionErrors) {
  EXPECT_DEATH({
    ЧебышёвSeries<double> p({}, t_min_, t_max_);
  }, "at least 0");
  EXPECT_DEATH({
    ЧебышёвSeries<double> p({1}, t_max_, t_min_);
  }, "not be empty");
}

TEST_F(ЧебышёвSeriesTest, EvaluationErrors) {
  EXPECT_DEATH({
    ЧебышёвSeries<double> p({1}, t_min_, t_max_);
    p.Evaluate(t_min_ - 10 * Second);
  }, ">= -1.1");
  EXPECT_DEATH({
    ЧебышёвSeries<double> p({1}, t_min_, t_max_);
    p.Evaluate(t_max_ + 10 * Second);
  }, "<= 1.1");
}

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

TEST_F(ЧебышёвSeriesDeathTest, SerializationError) {
  ЧебышёвSeries<Speed> v({1 * Metre / Second,
                          -2 * Metre / Second,
                          5 * Metre / Second},
                         t_min_, t_max_);
  ЧебышёвSeries<double> d({7, 8, -1}, t_min_, t_max_);

  EXPECT_DEATH({
    serialization::ЧебышёвSeries message;
    v.WriteToMessage(&message);
    ЧебышёвSeries<double>::ReadFromMessage(message);
  }, "has_double");
  EXPECT_DEATH({
    serialization::ЧебышёвSeries message;
    d.WriteToMessage(&message);
    ЧебышёвSeries<Speed>::ReadFromMessage(message);
  }, "has_quantity");
}

TEST_F(ЧебышёвSeriesTest, SerializationSuccess) {
  {
    serialization::ЧебышёвSeries message;
    ЧебышёвSeries<Speed> const v1({1 * Metre / Second,
                                   -2 * Metre / Second,
                                   5 * Metre / Second},
                                  t_min_, t_max_);
    v1.WriteToMessage(&message);
    EXPECT_EQ(3, message.coefficient_size());
    EXPECT_FALSE(message.coefficient(0).has_double_());
    EXPECT_TRUE(message.coefficient(0).has_quantity());
    EXPECT_EQ(0x7C01, message.coefficient(0).quantity().dimensions());
    EXPECT_EQ(1.0, message.coefficient(0).quantity().magnitude());
    EXPECT_TRUE(message.has_t_min());
    EXPECT_TRUE(message.t_min().has_scalar());
    EXPECT_TRUE(message.t_min().scalar().has_dimensions());
    EXPECT_TRUE(message.t_min().scalar().has_magnitude());
    EXPECT_EQ(-1.0, message.t_min().scalar().magnitude());
    EXPECT_TRUE(message.has_t_max());
    EXPECT_TRUE(message.t_max().has_scalar());
    EXPECT_TRUE(message.t_max().scalar().has_dimensions());
    EXPECT_TRUE(message.t_max().scalar().has_magnitude());
    EXPECT_EQ(3.0, message.t_max().scalar().magnitude());
    ЧебышёвSeries<Speed> const v2 =
        ЧебышёвSeries<Speed>::ReadFromMessage(message);
    EXPECT_EQ(v1, v2);
  }
  {
    serialization::ЧебышёвSeries message;
    ЧебышёвSeries<double> const d1({-1, 2, 5}, t_min_, t_max_);
    d1.WriteToMessage(&message);
    EXPECT_EQ(3, message.coefficient_size());
    EXPECT_TRUE(message.coefficient(0).has_double_());
    EXPECT_FALSE(message.coefficient(0).has_quantity());
    EXPECT_EQ(-1.0, message.coefficient(0).double_());
    EXPECT_TRUE(message.has_t_min());
    EXPECT_TRUE(message.t_min().has_scalar());
    EXPECT_TRUE(message.t_min().scalar().has_dimensions());
    EXPECT_TRUE(message.t_min().scalar().has_magnitude());
    EXPECT_EQ(-1.0, message.t_min().scalar().magnitude());
    EXPECT_TRUE(message.has_t_max());
    EXPECT_TRUE(message.t_max().has_scalar());
    EXPECT_TRUE(message.t_max().scalar().has_dimensions());
    EXPECT_TRUE(message.t_max().scalar().has_magnitude());
    EXPECT_EQ(3.0, message.t_max().scalar().magnitude());
    ЧебышёвSeries<double> const d2 =
        ЧебышёвSeries<double>::ReadFromMessage(message);
    EXPECT_EQ(d1, d2);
  }
}

}  // namespace numerics
}  // namespace principia
