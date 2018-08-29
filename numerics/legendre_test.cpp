
#include "numerics/legendre.hpp"

#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

class LegendreTest : public ::testing::Test {};

TEST_F(LegendreTest, P3P4) {
  auto const p3 = LegendrePolynomial<3, HornerEvaluator>();
  EXPECT_EQ(0, p3.Evaluate(0));
  EXPECT_EQ(1, p3.Evaluate(1));
  EXPECT_EQ(-1, p3.Evaluate(-1));

  auto const p4 = LegendrePolynomial<4, HornerEvaluator>();
  EXPECT_EQ(3.0 / 8.0, p4.Evaluate(0));
  EXPECT_EQ(1, p4.Evaluate(1));
  EXPECT_EQ(1, p4.Evaluate(-1));
  EXPECT_EQ(443.0 / 8.0, p4.Evaluate(2));
}

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia

