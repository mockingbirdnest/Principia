
#include "testing_utilities/statistics.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

using quantities::Inverse;
using quantities::Length;
using quantities::Speed;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;
using testing::Eq;

class StatisticsTest : public testing::Test {
 protected:
  StatisticsTest() {
    t_ = std::vector<Time>(population_size_);
    x_ = std::vector<Length>(population_size_);
    for (std::size_t i = 0; i < population_size_; ++i) {
      t_[i] = i / sampling_rate;
      x_[i] = t_[i] * v_ + x0_;
    }
  }

  std::size_t const population_size_ = 100;
  Inverse<Time> const sampling_rate = 8 / Second;
  Length const x0_ = - 12 * Metre;
  Speed const v_ = 42 * Metre / Second;
  std::vector<Time> t_;
  std::vector<Length> x_;
};

TEST_F(StatisticsTest, UniformPerfectlyCorrelated) {
  EXPECT_THAT(Mean(t_), Eq(((population_size_ - 1 + 0) / 2.0) / sampling_rate));
  EXPECT_THAT(Mean(x_), Eq((x_[population_size_ - 1] + x0_) / 2.0));
  EXPECT_THAT(Variance(t_),
              Eq(((population_size_ * population_size_ - 1) / 12.0) /
                 (sampling_rate * sampling_rate)));
  EXPECT_THAT(StandardDeviation(x_) * StandardDeviation(x_),
              AlmostEquals(Variance(x_), 1));
  EXPECT_THAT(StandardDeviation(x_) * StandardDeviation(t_),
              Eq(Covariance(x_, t_)));
  EXPECT_THAT(PearsonProductMomentCorrelationCoefficient(t_, x_), Eq(1));
  EXPECT_THAT(Slope(t_, x_), Eq(v_));
}

TEST_F(StatisticsTest, Uncorrelated) {
  std::vector<Time> t   = {0 * Second, 1 * Second, 1 * Second, 0 * Second};
  std::vector<Length> x = {0 * Metre, 0 * Metre, 1 * Metre, 1 * Metre};
  EXPECT_THAT(PearsonProductMomentCorrelationCoefficient(t, x), Eq(0));
}

TEST_F(StatisticsTest, NegativelyCorrelated) {
  std::vector<Time> t   = {0 * Second, 1 * Second, 1 * Second};
  std::vector<Length> x = {1 * Metre, 0 * Metre, 1 * Metre};
  EXPECT_THAT(PearsonProductMomentCorrelationCoefficient(t, x),
              AlmostEquals(-0.5, 1));
}

}  // namespace testing_utilities
}  // namespace principia
