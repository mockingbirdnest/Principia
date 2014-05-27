#include "testing_utilities/algebra.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/explicit_operators.hpp"

namespace principia {
namespace testing_utilities {

class AlgebraTest : public testing::Test {};

TEST_F(AlgebraTest, Group) {
  TestGroup(0, 42, -3, 2, Plus<int, int, int>, Minus<int, int>, 0);
}

}  // namespace testing_utilities
}  // namespace principia
