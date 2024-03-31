#include "base/cpuid.hpp"

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using ::testing::AnyOf;
using ::testing::HasSubstr;
using ::testing::Not;
using ::testing::Test;
using namespace principia::base::_cpuid;

class CPUIDTest : public Test {
 protected:
  CPUIDTest() {
    google::LogToStderr();
  }
};

TEST_F(CPUIDTest, Vendor) {
  // This mostly checks that we are getting something from CPUID, since it is
  // hard to expect things from the feature flags.  This could be expanded to an
  // AnyOf as needed if the tests are run on non-Intel processors.
  EXPECT_THAT(CPUVendorIdentificationString(),
              AnyOf("AuthenticAMD", "GenuineIntel"));
  LOG(INFO) << CPUVendorIdentificationString();
}

TEST_F(CPUIDTest, Brand) {
  EXPECT_THAT(ProcessorBrandString(),
              AnyOf(HasSubstr("AMD Ryzen"),
                    HasSubstr("Intel(R) Core(TM)"),
                    HasSubstr("Intel(R) Xeon(R)"),
                    HasSubstr("VirtualApple")));
  LOG(INFO) << ProcessorBrandString();
}

TEST_F(CPUIDTest, CPUFeatureFlags) {
  // We require Prescott or later.
  EXPECT_TRUE(CPUIDFeatureFlag::FPU.IsSet());
  EXPECT_TRUE(CPUIDFeatureFlag::SSE.IsSet());
  EXPECT_TRUE(CPUIDFeatureFlag::SSE2.IsSet());
  EXPECT_TRUE(CPUIDFeatureFlag::SSE3.IsSet());
  // Check that we don’t always return true.
  // We are not running these tests on a Pentium III, so we do not have the
  // Processor Serial Number feature.
  EXPECT_FALSE(CPUIDFeatureFlag::PSN.IsSet());
  EXPECT_THAT(
      CPUFeatures(),
      AllOf(HasSubstr("FPU"), HasSubstr("SSE2"), Not(HasSubstr("PSN"))));
  LOG(INFO) << CPUFeatures();
}

}  // namespace base
}  // namespace principia
