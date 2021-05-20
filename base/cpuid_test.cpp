#include "base/cpuid.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::Eq;
using ::testing::Test;

class CPUIDTest : public Test {};

TEST_F(CPUIDTest, Vendor) {
  // This mostly checks that we are getting something from CPUID, since it is
  // hard to expect things from the feature flags.  This could be expanded to an
  // AnyOf as needed if the tests are run on non-Intel processors.
  EXPECT_THAT(VendorIdentificationString(), Eq("GenuineIntel"));
}

TEST_F(CPUIDTest, FeatureFlags) {
  // We require Prescott or later.
  EXPECT_TRUE(HasCPUFeatures(FeatureFlags::FPU | FeatureFlags::SSE |
                             FeatureFlags::SSE2 | FeatureFlags::SSE3));
  // We develop on Sandy Bridge or later.
  EXPECT_TRUE(HasCPUFeatures(FeatureFlags::AVX));
  // Check that we don’t always return true.
  EXPECT_FALSE(HasCPUFeatures(FeatureFlags::NotUsed));
  EXPECT_FALSE(HasCPUFeatures(FeatureFlags::FPU | FeatureFlags::SSE |
                              FeatureFlags::SSE2 | FeatureFlags::SSE3 |
                              FeatureFlags::NotUsed));
}

}  // namespace base
}  // namespace principia
