
#include "quantities/wide.hpp"

#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

using si::Metre;

TEST(WideTest, Conversions) {
  Wide<Length> const w1(2 * Metre);
  __m128d m1 = ToM128D(w1);
  EXPECT_EQ(2, _mm_cvtsd_f64(m1));
  EXPECT_EQ(2, _mm_cvtsd_f64(_mm_unpackhi_pd(m1, m1)));

  Wide<double> const w2(3.14);
  __m128d m2 = ToM128D(w2);
  EXPECT_EQ(3.14, _mm_cvtsd_f64(m2));
  EXPECT_EQ(3.14, _mm_cvtsd_f64(_mm_unpackhi_pd(m2, m2)));

  Length const n3(3 * Metre);
  __m128d m3 = ToM128D(n3);
  EXPECT_EQ(3, _mm_cvtsd_f64(m3));
  EXPECT_EQ(3, _mm_cvtsd_f64(_mm_unpackhi_pd(m3, m3)));

  double const n4(2.71);
  __m128d m4 = ToM128D(n4);
  EXPECT_EQ(2.71, _mm_cvtsd_f64(m4));
  EXPECT_EQ(2.71, _mm_cvtsd_f64(_mm_unpackhi_pd(m4, m4)));
}

}  // namespace quantities
}  // namespace principia
