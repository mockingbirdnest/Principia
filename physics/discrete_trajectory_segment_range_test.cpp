#include "physics/discrete_trajectory_segment_range.hpp"

#include <vector>

#include "gtest/gtest.h"

namespace principia {
namespace physics {

TEST(DiscreteTrajectorySegmentRangeTest, Basic) {
  std::vector<int> const primes{2, 3, 5, 7, 11, 13, 17, 19, 23};
  DiscreteTrajectorySegmentRange<std::vector<int>::const_iterator> range(
      primes.begin(), primes.end());

  EXPECT_EQ(2, range.front());
  EXPECT_EQ(23, range.back());

  int sum = 0;
  int product = 1;
  for (int const prime : range) {
    sum += prime;
    product *= prime;
  }
  EXPECT_EQ(100, sum);
  EXPECT_EQ(223'092'870, product);

  EXPECT_EQ(9, primes.size());
  EXPECT_FALSE(primes.empty());
}

}  // namespace physics
}  // namespace principia
