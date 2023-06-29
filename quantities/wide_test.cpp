#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace quantities {

using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_matchers;

TEST(WideTest, Conversions) {
  __m128d m1 = ToM128D(2 * Metre);
  EXPECT_THAT(m1, SSEHighHalfIs(2));
  EXPECT_THAT(m1, SSELowHalfIs(2));

  __m128d m2 = ToM128D(3.14);
  EXPECT_THAT(m2, SSEHighHalfIs(3.14));
  EXPECT_THAT(m2, SSELowHalfIs(3.14));
}

}  // namespace quantities
}  // namespace principia
