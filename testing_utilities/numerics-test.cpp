#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/dimensionless.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

namespace {
struct World;
}  // namespace

class NumericsTest : public testing::Test {};

TEST_F(NumericsTest, ULPs) {
}

}  // namespace testing_utilities
}  // namespace principia
