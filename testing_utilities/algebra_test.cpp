#include "testing_utilities/algebra.hpp"

#include <functional>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace testing_utilities {

class AlgebraTest : public testing::Test {};

TEST_F(AlgebraTest, Group) {
  TestGroup(0, 42, -3, 2, std::plus<>(), std::negate<>(), 0, 0);
  TestGroup<double>(1.0, 42.0, -3.0, 2.0, std::multiplies<>(),
                    [](double const& x) { return 1 / x; }, 0, 0);
}

}  // namespace testing_utilities
}  // namespace principia
