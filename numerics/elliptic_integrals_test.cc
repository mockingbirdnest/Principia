
#include <filesystem>
#include <vector>
#include <utility>

#include "numerics/elliptic_integrals.hpp"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"

namespace principia {

using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::ReadFromTabulatedData;

namespace numerics {

class EllipticIntegralsTest : public ::testing::Test {};

// TODO(phl): This is far from covering all the code paths.  In particular, it
// doesn't seem to go through Elcbdj.
TEST_F(EllipticIntegralsTest, Xelbdj) {
  auto const xeldbj_expected =
      ReadFromTabulatedData(SOLUTION_DIR / "numerics" / "xelbdj.proto.txt");
  constexpr double PI = 3.1415926535897932384626433;
  constexpr double PIHALF = 1.5707963267948966192313216916398;
  double dnc, dmc, dphi, nc, nn, mc, mm, phi, phic, b, d, j;
  int lend, kend, iend;

  lend = 5;
  kend = 5;
  iend = 4;
  dnc = 1.0 / static_cast<double>(lend - 1);
  dmc = 1.0 / static_cast<double>(kend - 1);
  dphi = PIHALF / static_cast<double>(iend);
  std::printf(
      "%10s%10s%10s%25s%25s%25s\n", "n", "m", "phi/PI", "elb", "eld", "elj");
  int expected_index = 0;
  for (int l = 1; l <= lend; ++l) {
    nc = static_cast<double>(l - 1) * dnc;
    if (nc <= 1.05e-8) {
      nc = 1.05e-8;
    }
    float rnc = static_cast<float>(nc);
    if (rnc <= 2.44e-4) {
      nc = 2.44e-4;
    }
    nn = 1.0 - nc;
    for (int k = 1; k <= kend; ++k) {
      std::printf("\n");
      mc = static_cast<double>(k - 1) * dmc;
      if (mc <= 0.0) {
        mc = 1.21e-32;
      }
      mm = 1.0 - mc;
      for (int i = 0; i <= iend; ++i) {
        phi = dphi * static_cast<double>(i);
        phic = dphi * static_cast<double>(iend - i);
        Elbdj(phi, phic, nn, mc, b, d, j);
        std::printf("%10.5f%10.5f%10.5f%25.15f%25.15f%25.15f\n",
                    nn,
                    mm,
                    phi / PI,
                    b,
                    d,
                    j);

        auto const& expected_entry = xeldbj_expected.entry(expected_index);
        auto const expected_argument_n = expected_entry.argument(0);
        auto const expected_argument_m = expected_entry.argument(1);
        auto const expected_argument_phi_over_pi = expected_entry.argument(2);
        auto const expected_value_b = expected_entry.value(0);
        auto const expected_value_d = expected_entry.value(1);
        auto const expected_value_j = expected_entry.value(2);
        EXPECT_THAT(nn, IsNear(expected_argument_n, 1.001));
        EXPECT_THAT(mm, IsNear(expected_argument_m, 1.001));
        EXPECT_THAT(phi / PI, IsNear(expected_argument_phi_over_pi, 1.001));
        EXPECT_THAT(b, AlmostEquals(expected_value_b, 0, 8));
        EXPECT_THAT(d, AlmostEquals(expected_value_d, 0, 97));
        EXPECT_THAT(j, AlmostEquals(expected_value_j, 0, 135));
        ++expected_index;
      }
    }
  }
}

}  // namespace numerics
}  // namespace principia
