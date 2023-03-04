#include "testing_utilities/algebra.hpp"

#include <functional>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define PRINCIPIA_COMPILER_MSVC_HANDLES_TEST_GROUP \
  (!PRINCIPIA_COMPILER_MSVC || !(_MSC_FULL_VER == 193'431'942))

namespace principia {
namespace testing_utilities {
#if PRINCIPIA_COMPILER_MSVC_HANDLES_TEST_GROUP
using namespace principia::testing_utilities::_algebra;
#else
namespace _algebra {
namespace internal {
#endif

class AlgebraTest : public testing::Test {};

TEST_F(AlgebraTest, Group) {
  TestGroup(0, 42, -3, 2, std::plus<>(), std::negate<>(), 0, 0);
  TestGroup<double>(1.0, 42.0, -3.0, 2.0, std::multiplies<>(),
                    [](double const& x) { return 1 / x; }, 0, 0);
}

#if !PRINCIPIA_COMPILER_MSVC_HANDLES_TEST_GROUP
}  // namespace internal
}  // namespace _algebra
#endif
}  // namespace testing_utilities
}  // namespace principia
