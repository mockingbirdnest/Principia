
#include "numerics/apodization.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using geometry::Instant;
using quantities::si::Second;
using testing_utilities::AlmostEquals;

class ApodizationTest : public ::testing::Test {
 protected:
  ApodizationTest() : t1_(t0_ - 1 * Second), t2_(t0_ + 2 * Second) {}

  Instant const t0_;
  Instant const t1_;
  Instant const t2_;
};

TEST_F(ApodizationTest, Dirichlet) {
  auto a = apodization::Dirichlet<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(1, AlmostEquals(a.Evaluate(t0_), 0));
  EXPECT_THAT(1, AlmostEquals(a.Evaluate(t0_ + 0.5 * Second), 0));
  EXPECT_THAT(1, AlmostEquals(a.Evaluate(t0_ + 1 * Second), 0));
}

}  // namespace numerics
}  // namespace principia
