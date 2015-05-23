
#include "numerics/чебышёв_series.hpp"

#include <algorithm>
#include <cmath>

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using geometry::Instant;
using quantities::Length;
using quantities::Speed;
using si::Metre;
using si::Second;
using testing_utilities::AbsoluteError;
using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Gt;
using ::testing::Lt;

namespace numerics {

class ЧебышёвSeriesTest : public ::testing::Test {
 protected:
  ЧебышёвSeriesTest()
      : t_min_(-1 * Second),
        t_max_(3 * Second) {}

  // TODO(phl): Support for evaluating the derivative.
  void NewhallApproximationErrors(
      std::function<Length(Instant const)> length_function,
      std::function<Speed(Instant const)> speed_function,
      not_null<std::vector<Length>*> const length_absolute_errors) {
    std::vector<Length> lengths;
    std::vector<Speed> speeds;
    for (Instant t = t_min_; t <= t_max_; t += 0.5 * Second) {
      lengths.push_back(length_function(t));
      speeds.push_back(speed_function(t));
    }

    length_absolute_errors->clear();
    for (int degree = 3; degree <= 17; ++degree) {
      ЧебышёвSeries<Length> const approximation =
          ЧебышёвSeries<Length>::NewhallApproximation(
              degree, lengths, speeds, t_min_, t_max_);
      Length absolute_error;
      for (Instant t = t_min_; t <= t_max_; t += 0.05 * Second) {
        Length const expected_length = length_function(t);
        Length const actual_length = approximation.Evaluate(t);
        absolute_error = std::max(absolute_error,
                                  AbsoluteError(expected_length, actual_length));
      }
      if (!length_absolute_errors->empty()) {
        //CHECK_LE(absolute_error, length_absolute_errors->back());
      }
      length_absolute_errors->push_back(absolute_error);
    }
  }

  // A helper that splits an array in two chunks and applies distinct matchers
  // to the chunks.  Necessary because ElementsAre only supports 10 elements and
  // ElementsAreArray does not support matches as arguments.
  template<typename Matcher1, typename Matcher2>
  void ExpectMultipart(std::vector<Length> const& v,
                       Matcher1 const& matcher1,
                       Matcher2 const& matcher2) {
    std::vector<Length> v_0_9(10);
    std::copy(v.begin(), v.begin() + 10, v_0_9.begin());
    EXPECT_THAT(v_0_9, matcher1);
    std::vector<Length> v_10_end(v.size() - 10);
    std::copy(v.begin() + 10, v.end(), v_10_end.begin());
    EXPECT_THAT(v_10_end, matcher2);
  }

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

TEST_F(ЧебышёвSeriesTest, T2Dimension) {
  ЧебышёвSeries<Length> t2({0 * Metre, 0 * Metre, 1 * Metre}, t_min_, t_max_);
  EXPECT_EQ(1 * Metre, t2.Evaluate(Instant(-1 * Second)));
  EXPECT_EQ(-1 * Metre, t2.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1 * Metre, t2.Evaluate(Instant(3 * Second)));
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

TEST_F(ЧебышёвSeriesTest, NewhallApproximation) {
  std::vector<Length> length_absolute_errors;

  auto near = [](Length const& length) {
    return AllOf(Gt(0.9 * length), Lt(length));
  };

  {
    auto length_function = [this](Instant const t) -> Length {
      return 2 * Metre *
             std::sin((t - t_min_) / (0.3 * Second)) *
             std::exp((t - t_min_) / (1 * Second));
    };
    auto speed_function = [this](Instant const t) -> Speed {
      return ((2 * Metre) / (0.3 * Second) *
                   std::cos((t - t_min_) / (0.3 * Second)) +
              (2 * Metre / Second) *
                   std::sin((t - t_min_) / (0.3 * Second))) *
                       std::exp((t - t_min_) / (1 * Second));
    };

    NewhallApproximationErrors(length_function,
                               speed_function,
                               &length_absolute_errors);
  }

  ExpectMultipart(length_absolute_errors,
                  ElementsAre(near(1.7E2 * Metre),
                              near(4.7E1 * Metre),
                              near(4.3E1 * Metre),
                              near(3.8E1 * Metre),
                              near(1.5E1 * Metre),
                              near(6.3 * Metre),
                              near(4.9 * Metre),
                              near(6.5E-1 * Metre),
                              near(2.0E-1 * Metre),
                              near(7.9E-2 * Metre)),
                  ElementsAre(near(1.3E-2 * Metre),
                              near(1.6E-2 * Metre),
                              near(4.3E-3 * Metre),
                              near(1.7E-3 * Metre),
                              near(7.6E-4 * Metre)));

  {
    auto length_function = [this](Instant const t) -> Length {
      return 5 * Metre * ((t - t_min_) / (0.3 * Second) +
                          std::pow((t - t_min_) / (4 * Second), 7));
    };
    auto speed_function = [this](Instant const t) -> Speed {
      return 5 * Metre * (1 / (0.3 * Second) +
                          (7 / (4 * Second)) *
                              std::pow((t - t_min_) / (4 * Second), 6));
    };

    NewhallApproximationErrors(length_function,
                               speed_function,
                               &length_absolute_errors);
  }

  ExpectMultipart(length_absolute_errors,
                  ElementsAre(near(2.0 * Metre),
                              near(2.9E-1 * Metre),
                              near(3.6E-2 * Metre),
                              near(2.3E-3 * Metre),
                              near(2.9E-14 * Metre),
                              near(2.9E-14 * Metre),
                              near(2.9E-14 * Metre),
                              near(2.0E-14 * Metre),
                              near(2.4E-14 * Metre),
                              near(2.9E-14 * Metre)),
                  ElementsAre(near(2.2E-14 * Metre),
                              near(1.5E-14 * Metre),
                              near(2.2E-14 * Metre),
                              near(1.6E-14 * Metre),
                              near(4.3E-14 * Metre)));
}

}  // namespace numerics
}  // namespace principia
