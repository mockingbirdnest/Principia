#include "numerics/lattices.hpp"

#include "boost/multiprecision/cpp_int.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {
namespace _lattices {

using namespace boost::multiprecision;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_lattices;
using namespace principia::testing_utilities::_almost_equals;

class LatticesTest : public ::testing::Test {};

// [HPS14], example 7.75.
TEST_F(LatticesTest, LLL_Example_7_75) {
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

TEST_F(LatticesTest, LLL_Rational) {
  FixedMatrix<cpp_rational, 5, 4> l({
      45,         0, 214695880217044191, 401754430875619365,
       0, 188743680,          187081485,           -6248472,
       0,         0,                  0,                  0,
       0,         0,                  3,                  0,
       0,         0,                  0,                  3});

  auto const reduced = LenstraLenstraLovász(l);
  EXPECT_EQ(reduced,
            (FixedMatrix<cpp_rational, 5, 4>({ 45,    6,   -18,   15,
                                                0,   45, -1200,  348,
                                                0,    0,     0,    0,
                                                0, 1083,   336, -660,
                                                0,  165,  -180, 1263})));
}

TEST_F(LatticesTest, NS_Example_7_75) {
  FixedMatrix<double, 6, 6> l({19, 15, 43, 20,  0, 48,
                                2, 42, 15, 44, 48, 33,
                               32, 11,  0, 44, 35, 32,
                               46,  0, 24,  0, 16,  9,
                                3,  3,  4, 18, 31,  1,
                               33, 24, 16, 15, 31, 29});

  auto const reduced = NguyễnStehlé(l);
  // Note that, even taking into account permutations and sign changes, this
  // reduction is not the same as the one obtained by LLL: the last vector
  // differs.  It appears to be a little longer and somewhat more orthogonal
  // than the LLL one.
  EXPECT_THAT(reduced,
              AlmostEquals(
                  (FixedMatrix<double, 6, 6>({ -7, -20,  6,  5,  -7,  -9,
                                               12,   4,  7,  2,  -4, -19,
                                                8,  -9, 20, 33,   9,   8,
                                               -4,  16, 21,  0,  11,   6,
                                              -19,  13, -8, 15,  -1, -29,
                                               -9,  16, 12, -9, -31,  10})),
                  0));
}

TEST_F(LatticesTest, NS_Int) {
  FixedMatrix<cpp_int, 5, 4> l({
      45,         0, 214695880217044191, 401754430875619365,
       0, 188743680,          187081485,           -6248472,
       0,         0,                  0,                  0,
       0,         0,                  3,                  0,
       0,         0,                  0,                  3});

  auto const reduced = NguyễnStehlé(l);
  EXPECT_EQ(reduced,
            (FixedMatrix<cpp_int, 5, 4>({ 45,    6,   -18,   15,
                                           0,   45, -1200,  348,
                                           0,    0,     0,    0,
                                           0, 1083,   336, -660,
                                           0,  165,  -180, 1263})));
}

TEST_F(LatticesTest, NoProgress) {
  FixedMatrix<cpp_int, 5, 4> l(
      { 212,  74278, -24016,  cpp_int("2236050964407517784120"),
       -116, 135836, -23372, cpp_int("-1223499584298453127159"),
          0,      0,      0,                                  0,
          3,  -3513, 784686,    cpp_int("31642230628408270524"),
          0,      0,      0,                                  3});
  auto const reduced = NguyễnStehlé(l);
  EXPECT_EQ(reduced,
            (FixedMatrix<cpp_int, 5, 4>({ 0,  212,  72582,  26440,
                                          1, -113, 133330,  52572,
                                          0,    0,      0,      0,
                                         -6,  -15,  17067, 164088,
                                          3,    9, -10302, 310656})));
}

}  // namespace _lattices
}  // namespace numerics
}  // namespace principia
