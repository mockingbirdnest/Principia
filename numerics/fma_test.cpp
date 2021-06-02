
#include "numerics/fma.hpp"

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using base::CPUFeatureFlags;
using base::HasCPUFeatures;
using testing_utilities::AlmostEquals;

class FMATest : public testing::Test {
 protected:
  void SetUp() override {
    // Note that we test even if |UseHardwareFMA| is false, i.e., even in debug.
    if (!CanEmitFMAInstructions || !HasCPUFeatures(CPUFeatureFlags::FMA)) {
      GTEST_SKIP() << "Cannot test FMA on a machine without FMA";
    }
  }
};

TEST_F(FMATest, FMA) {
  EXPECT_THAT(FusedMultiplyAdd(0.2, 5.1, 1.2),
              AlmostEquals(0.2 * 5.1 + 1.2, 0));
  EXPECT_THAT(FusedMultiplySubtract(0.2, 5.1, 1.2),
              AlmostEquals(0.2 * 5.1 - 1.2, 1));
  EXPECT_THAT(FusedNegatedMultiplyAdd(0.2, 5.1, 1.2),
              AlmostEquals(-0.2 * 5.1 + 1.2, 1));
  EXPECT_THAT(FusedNegatedMultiplySubtract(0.2, 5.1, 1.2),
              AlmostEquals(-0.2 * 5.1 - 1.2, 0));
}

}  // namespace numerics
}  // namespace principia
