
#include "numerics/legendre.hpp"

#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

class LegendreTest : public ::testing::Test {
 protected:
  LegendrePolynomial<2, HornerEvaluator> p2_;
};

TEST_F(LegendreTest, Smoke) {}

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia

