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
using testing::Eq;

namespace principia {
namespace testing_utilities {

class VanishesBeforeTest : public testing::Test {};

TEST_F(VanishesBeforeTest, R3Element) {
  R3Element<double> r1({1.0 + 1.0E-12, 1.0E-10, 3.5});
  EXPECT_THAT(r1, Componentwise(AlmostEquals(1.0, 45044),
                                VanishesBefore(1.0, 450360),
                                Eq(3.5)));
}

}  // namespace testing_utilities
}  // namespace principia
