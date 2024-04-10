#include "functions/accurate_table_generator.hpp"

#include <iomanip>

#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {

using ::testing::AnyOf;
using namespace boost::multiprecision;
using namespace principia::functions::_accurate_table_generator;
using namespace principia::functions::_multiprecision;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;

class AccurateTableGeneratorTest : public ::testing::Test {};

TEST_F(AccurateTableGeneratorTest, Sin5) {
  auto const x = ExhaustiveSearch<5>({Sin}, 5.0 / 128.0);
  EXPECT_EQ(x, cpp_rational(2814749767106647, 72057594037927936));
  EXPECT_THAT(static_cast<double>(x),
              RelativeErrorFrom(5.0 / 128.0, IsNear(3.1e-14_(1))));

  // The stupid language doesn't allow printing a float in binary.  So to verify
  // that the |Sin| has zeroes in the right place we first produce a hex string
  // with significantly more bits than the mantissa.  Then we check that:
  // 1. The leading digit has its leading bit set to 1 (for alignment).
  // 2. There are 12 digits that we ignore after it for a total of 52 bits.
  // 3. The next two nibbles are 16#83# or 2#1000_0011#, where the leading 1 is
  //    part of the mantissa, and the following 5 zeroes are the ones we are
  //    after.
  // This validates that we actually shoot zeroes at the right place.
  std::string const s =
      (std::stringstream() << std::hex
                           << static_cast<cpp_int>(floor(ldexp(Sin(x), 64))))
          .str();
  EXPECT_THAT(s[0], AnyOf('8', '9', 'a', 'b', 'c', 'd', 'e', 'f'));
  EXPECT_EQ(s[13], '8');
  EXPECT_EQ(s[14], '3');
}

}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
