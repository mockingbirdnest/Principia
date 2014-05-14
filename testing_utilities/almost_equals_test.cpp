#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "geometry/grassmann.hpp"
#include "quantities/dimensionless.hpp"
#include "quantities/quantities.hpp"
#include "quantities/uk.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

using geometry::R3Element;
using geometry::Vector;
using quantities::Dimensionless;
using quantities::Length;
using testing::Eq;
using testing::Not;
using uk::Foot;

class AlmostEqualsTest : public testing::Test {
 protected:
};

TEST_F(AlmostEqualsTest, AlmostButNotQuiteEquals) {
  EXPECT_THAT(Dimensionless(1), AlmostEquals(1));
  EXPECT_THAT(Dimensionless(1.01), Not(AlmostEquals(1)));
  Dimensionless not_quite_one = 0;
  for (int i = 1; i <= 20; ++i) {
    not_quite_one += 0.05;
  }
  EXPECT_THAT(not_quite_one, Not(Eq(1)));
  EXPECT_THAT(not_quite_one, AlmostEquals(1));
  EXPECT_THAT(not_quite_one - 1, Not(AlmostEquals(0)));
  EXPECT_THAT(not_quite_one - 1, AlmostVanishesBefore(1));
  EXPECT_THAT(not_quite_one - 1, Not(AlmostVanishesBefore(.2)));
}

TEST_F(AlmostEqualsTest, ApproximationMatcher) {
  EXPECT_THAT(Dimensionless(2.19), Approximates(2, 0.1));
  EXPECT_THAT(Dimensionless(2.21), Not(Approximates(2, 0.1)));
}

TEST_F(AlmostEqualsTest, DimensionfulAlmostButNotQuiteEquals) {
  EXPECT_THAT(1 * Foot, AlmostEquals(1 * Foot));
  EXPECT_THAT(1.01 * Foot, Not(AlmostEquals(1 * Foot)));
  Length not_quite_one = 0 * Foot;
  for (int i = 1; i <= 20; ++i) {
    not_quite_one += 0.05 * Foot;
  }
  EXPECT_THAT(not_quite_one, Not(Eq(1 * Foot)));
  EXPECT_THAT(not_quite_one, AlmostEquals(1 * Foot));
  EXPECT_THAT(not_quite_one - 1 * Foot, Not(AlmostEquals(0 * Foot)));
  EXPECT_THAT(not_quite_one - 1 * Foot, AlmostVanishesBefore(1 * Foot));
  EXPECT_THAT(not_quite_one - 1 * Foot, Not(AlmostVanishesBefore(.2 * Foot)));
}

TEST_F(AlmostEqualsTest, DimensionfulApproximationMatcher) {
  EXPECT_THAT(2.19 * Foot, Approximates(2 * Foot, 0.1));
  EXPECT_THAT(2.21 * Foot, Not(Approximates(2 * Foot, 0.1)));
}

struct World;
/*
TEST_F(AlmostEqualsTest, Vectors) {
  auto v1 = Vector<Dimensionless, World>(R3Element<Dimensionless>(1, 2, 3));
  auto v2 = v1;
  EXPECT_THAT(v1, AlmostEquals(v2));
}
*/

}  // namespace testing_utilities
}  // namespace principia
