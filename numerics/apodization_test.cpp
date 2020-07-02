
#include "numerics/apodization.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using geometry::Instant;
using quantities::si::Second;

class ApodizationTest : public ::testing::Test {
 protected:
   ApodizationTest() : t1_(t0_ - 1 * Second), t2_(t0_ + 2 * Second) {}

  Instant const t0_;
  Instant const t1_;
  Instant const t2_;
};

TEST_F(ApodizationTest, Dirichlet) {
  auto a = apodization::Dirichlet<HornerEvaluator>(t1_, t2_);
  EXPECT_EQ(1, a.Evaluate(t0_ + 0.5 * Second));
}

}  // namespace numerics
}  // namespace principia
