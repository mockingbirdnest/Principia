#include "numerics/lattices.hpp"

#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {
namespace _lattices {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_lattices;
using namespace principia::testing_utilities::_almost_equals;

class LatticesTest : public ::testing::Test {};

// [HPS14], example 7.75.
TEST_F(LatticesTest, Example_7_75) {
  FixedMatrix<double, 6, 6> l({19, 15, 43, 20,  0, 48,
                                2, 42, 15, 44, 48, 33,
                               32, 11,  0, 44, 35, 32,
                               46,  0, 24,  0, 16,  9,
                                3,  3,  4, 18, 31,  1,
                               33, 24, 16, 15, 31, 29});

  auto const reduced = LenstraLenstraLovász(l);
  EXPECT_THAT(reduced,
              AlmostEquals(
                  (FixedMatrix<double, 6, 6>({  7, -20,  5,  -6, -10,   7,
                                              -12,   4,  2,  -7, -24,   4,
                                               -8,  -9, 33, -20,  21,  -9,
                                                4,  16,  0, -21, -15, -11,
                                               19,  13, 15,   8,  -6,   1,
                                                9,  16, -9, -12, -11,  31})),
                  0));
}

}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
