
#include "quantities/wide.hpp"

#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

using si::Metre;

TEST(WideTest, ConstructionConversion) {
  Wide<Length> const w1(2 * Metre);
  __m128d m1 = w1.m128d();
  EXPECT_EQ(2, m1.m128d_f64[0]);
  EXPECT_EQ(2, m1.m128d_f64[1]);

  Wide<double> const w2(3.14);
  __m128d m2 = w2.m128d();
  EXPECT_EQ(3.14, m2.m128d_f64[0]);
  EXPECT_EQ(3.14, m2.m128d_f64[1]);
}

}  // namespace quantities
}  // namespace principia
