
#include "numerics/apodization.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace numerics {

using geometry::Instant;
using quantities::Sqrt;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;

class ApodizationTest : public ::testing::Test {
 protected:
  ApodizationTest()
      : t1_(t0_ - 1 * Second),
        t2_(t0_ + 2 * Second),
        mid_(t0_ + 0.5 * Second) {}

  Instant const t0_;
  Instant const t1_;
  Instant const t2_;
  Instant const mid_;
};

TEST_F(ApodizationTest, Dirichlet) {
  auto a = apodization::Dirichlet<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a.Evaluate(t1_), AlmostEquals(1, 0));
  EXPECT_THAT(a.Evaluate(t0_), AlmostEquals(1, 0));
  EXPECT_THAT(a.Evaluate(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a.Evaluate(t2_), AlmostEquals(1, 0));
}

TEST_F(ApodizationTest, Sine) {
  auto a = apodization::Sine<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a.Evaluate(t1_), AlmostEquals(0, 0));
  EXPECT_THAT(a.Evaluate(t0_), AlmostEquals(Sqrt(3) / 2, 0));
  EXPECT_THAT(a.Evaluate(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a.Evaluate(t2_), VanishesBefore(1, 1));
}

TEST_F(ApodizationTest, Hahn) {
  auto a = apodization::Hahn<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a.Evaluate(t1_), AlmostEquals(0, 0));
  EXPECT_THAT(a.Evaluate(t0_), AlmostEquals(0.75, 1));
  EXPECT_THAT(a.Evaluate(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a.Evaluate(t2_), VanishesBefore(1, 0));
}

}  // namespace numerics
}  // namespace principia
