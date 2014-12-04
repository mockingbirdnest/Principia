#include "testing_utilities/componentwise.hpp"

#include <limits>

#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

using principia::geometry::R3Element;
using testing::AllOf;

namespace principia {
namespace testing_utilities {

class VanishesBeforeTest : public testing::Test {};

TEST_F(VanishesBeforeTest, R3Element) {
  R3Element<double> r1({1.0 + 1.0E-12, 1.0E-10, 1.0E-11});
  EXPECT_THAT(r1, Componentwise(AlmostEquals(1.0, 4),
                                VanishesBefore(1.0, 2),
                                VanishesBefore(2.0, 1)));
}

}  // namespace testing_utilities
}  // namespace principia
