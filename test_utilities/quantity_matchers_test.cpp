#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "quantities/dimensionless.hpp"
#include "test_utilities/quantity_matchers.hpp"

namespace principia {
namespace test_utilities {

using quantities::Dimensionless;

class QuantityMatchersTest : public testing::Test {
 protected:
};

TEST_F(QuantityMatchersTest, Foo) {
  EXPECT_THAT(Dimensionless(1), AlmostEquals(Dimensionless(1)));
  EXPECT_THAT(Dimensionless(2), AlmostEquals(Dimensionless(1)));
}

}  // namespace test_utilities
}  // namespace principia
