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
  EXPECT_THAT(CPUVendorIdentificationString(), Eq("GenuineIntel"));
}

TEST_F(CPUIDTest, CPUFeatureFlags) {
  // We require Prescott or later.
  EXPECT_TRUE(HasCPUFeatures(CPUFeatureFlags::FPU | CPUFeatureFlags::SSE |
                             CPUFeatureFlags::SSE2 | CPUFeatureFlags::SSE3));
  // Check that we don’t always return true.
  // We are not running these tests on a Pentium III, so we do not have the
  // Processor Serial Number feature.
  EXPECT_FALSE(HasCPUFeatures(CPUFeatureFlags::PSN));
  EXPECT_FALSE(HasCPUFeatures(CPUFeatureFlags::FPU | CPUFeatureFlags::SSE |
                              CPUFeatureFlags::SSE2 | CPUFeatureFlags::SSE3 |
                              CPUFeatureFlags::PSN));
}

}  // namespace base
}  // namespace principia
