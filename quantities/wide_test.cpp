
#include "quantities/wide.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace quantities {

using si::Metre;
using testing_utilities::SSEHighHalfIs;
using testing_utilities::SSELowHalfIs;
using ::testing::Eq;

TEST(WideTest, Conversions) {
  Wide<Length> const w1(2 * Metre);
  __m128d m1 = ToM128D(w1);
  EXPECT_THAT(m1, SSEHighHalfIs(2));
  EXPECT_THAT(m1, SSELowHalfIs(2));

  Wide<double> const w2(3.14);
  __m128d m2 = ToM128D(w2);
  EXPECT_THAT(m2, SSEHighHalfIs(3.14));
  EXPECT_THAT(m2, SSELowHalfIs(3.14));

  Length const n3(3 * Metre);
  __m128d m3 = ToM128D(n3);
  EXPECT_THAT(m3, SSEHighHalfIs(3));
  EXPECT_THAT(m3, SSELowHalfIs(3));

  double const n4(2.71);
  __m128d m4 = ToM128D(n4);
  EXPECT_THAT(m4, SSEHighHalfIs(2.71));
  EXPECT_THAT(m4, SSELowHalfIs(2.71));
}

}  // namespace quantities
}  // namespace principia
