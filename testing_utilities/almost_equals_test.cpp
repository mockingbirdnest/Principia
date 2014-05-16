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
using testing::Ne;
using testing::Eq;
using testing::Not;
using uk::Foot;

class AlmostEqualsTest : public testing::Test {
 protected:
};

struct World;
TEST_F(AlmostEqualsTest, Vectors) {
  auto const v1 = Vector<Length, World>(
      R3Element<Length>(1 * Foot, 2 * Foot, 3 * Foot));
  Vector<Length, World> const v2 = v1;
  EXPECT_THAT(v2, AlmostEquals(v1));
  EXPECT_THAT(2 * v2, Not(AlmostEquals(v1)));
  Vector<Length, World> const δv = v1 / 100;
  auto v_accumulated = Vector<Length, World>(R3Element<Length>());
  for (int i = 1; i <= 100; ++i) {
    v_accumulated += δv;
  }
  EXPECT_THAT(v_accumulated, Ne(v1));
  EXPECT_THAT(v_accumulated, AlmostEquals(v1, 14));
}

}  // namespace testing_utilities
}  // namespace principia
