#include "testing_utilities/algebra.hpp"

#include <functional>

#include "base/algebra.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace testing_utilities {
namespace _algebra {
namespace internal {

using namespace base::_algebra;

class AlgebraTest : public testing::Test {};

static_assert(std::equality_comparable<IntegerModulo<2>>);
static_assert(!std::totally_ordered<IntegerModulo<2>>);
static_assert(additive_group<IntegerModulo<2>>);
static_assert(ring<IntegerModulo<2>>);
static_assert(field<IntegerModulo<2>>);
static_assert(!real_vector_space<IntegerModulo<2>>);

static_assert(std::equality_comparable<IntegerModulo<4>>);
static_assert(!std::totally_ordered<IntegerModulo<4>>);
static_assert(additive_group<IntegerModulo<4>>);
static_assert(ring<IntegerModulo<4>>);
static_assert(!field<IntegerModulo<4>>);

TEST_F(AlgebraTest, Group) {
  TestGroup(0, 42, -3, 2, std::plus<>(), std::negate<>(), 0, 0);
  TestGroup<double>(1.0, 42.0, -3.0, 2.0, std::multiplies<>(),
                    [](double const& x) { return 1 / x; }, 0, 0);
}

}  // namespace internal
}  // namespace _algebra
}  // namespace testing_utilities
}  // namespace principia
