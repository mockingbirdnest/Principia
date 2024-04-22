#include "numerics/lattices.hpp"

#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace numerics {
namespace _lattices {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_lattices;

class LatticesTest : public ::testing::Test {};

TEST_F(LatticesTest, Example_7_75) {
  FixedMatrix<double, 6, 6> l({19, 15, 43, 20,  0, 48,
                                2, 42, 15, 44, 48, 33,
                                32, 11,  0, 44, 35, 32,
                                46,  0, 24,  0, 16,  9,
                                3,  3,  4, 18, 31,  1,
                                33, 24, 16, 15, 31, 29});
  auto const reduced = LenstraLenstraLovász(l);
  LOG(ERROR)<<reduced;
}

}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
