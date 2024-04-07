#include <cmath>
#include <random>

#include "functions/multiprecision.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {

using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;
using namespace boost::multiprecision;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_numerics;

class StdAccuracyTest : public ::testing::Test {};

TEST_F(StdAccuracyTest, Sin) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> angle_distribution(0, π / 4);
  for (int i = 0; i < 1000; ++i) {
    double α = angle_distribution(random);
    EXPECT_THAT(RelativeError(std::sin(α),
                              static_cast<double>(Sin(cpp_rational(α)))),
                AllOf(Ge(0.0), Le(2.2e-16)));
  }
}

}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
