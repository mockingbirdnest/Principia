#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace integrators {

namespace {

struct SimpleHarmonicMotionTestInstance {};

std::vector<SimpleHarmonicMotionTestInstance> Instances() {}

}  // namespace

class SymmetricLinearMultistepIntegratorTest
    : public ::testing::TestWithParam<SimpleHarmonicMotionTestInstance> {};

}  // namespace integrators
}  // namespace principia
