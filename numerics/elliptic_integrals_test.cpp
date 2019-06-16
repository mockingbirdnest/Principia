
#include <filesystem>
#include <vector>
#include <utility>

#include "numerics/elliptic_integrals.hpp"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"

#include "quantities/quantities.hpp"

namespace principia {

using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::ReadFromTabulatedData;

namespace numerics {

class EllipticIntegralsTest : public ::testing::Test {};

// The test values found in Fukushima's xelbdj_all.txt file.
// NOTE(phl): This is far from covering all the code paths.  In particular, it
// doesn't seem to go through Elcbdj.
TEST_F(EllipticIntegralsTest, Xelbdj) {
  auto const xeldbj_expected =
      ReadFromTabulatedData(SOLUTION_DIR / "numerics" / "xelbdj.proto.txt");
  double Δnc, Δmc, Δφ, nc, nn, mc, mm, φ, b, d, j;
  int lend, kend, iend;

  lend = 5;
  kend = 5;
  iend = 4;
  Δnc = 1.0 / static_cast<double>(lend - 1);
  Δmc = 1.0 / static_cast<double>(kend - 1);
  Δφ = (π / 2) / static_cast<double>(iend);
  std::printf(
      "%10s%10s%10s%25s%25s%25s\n", "n", "m", "φ / π", "elb", "eld", "elj");
  int expected_index = 0;
  for (int l = 1; l <= lend; ++l) {
    nc = static_cast<double>(l - 1) * Δnc;
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
      mc = static_cast<double>(k - 1) * Δmc;
      if (mc <= 0.0) {
        mc = 1.21e-32;
      }
      mm = 1.0 - mc;
      for (int i = 0; i <= iend; ++i) {
        φ = Δφ * static_cast<double>(i);
        FukushimaEllipticBDJ(φ, nn, mc, b, d, j);
        std::printf("%10.5f%10.5f%10.5f%25.15f%25.15f%25.15f\n",
                    nn,
                    mm,
                    φ / π,
                    b,
                    d,
                    j);

        auto const& expected_entry = xeldbj_expected.entry(expected_index);
        auto const expected_argument_n = expected_entry.argument(0);
        auto const expected_argument_m = expected_entry.argument(1);
        auto const expected_argument_φ_over_π = expected_entry.argument(2);
        auto const expected_value_b = expected_entry.value(0);
        auto const expected_value_d = expected_entry.value(1);
        auto const expected_value_j = expected_entry.value(2);
        EXPECT_THAT(nn, IsNear(expected_argument_n, 1.001));
        EXPECT_THAT(mm, IsNear(expected_argument_m, 1.001));
        EXPECT_THAT(φ / π, IsNear(expected_argument_φ_over_π, 1.001));
        EXPECT_THAT(b, AlmostEquals(expected_value_b, 0, 8));

        // TODO(phl): xelbdj_all.txt enshrines values that are incorrect for
        // m close to 1 and φ close to (the machine representation of) π / 2.
        // That seems to stem from the computation of the complementary angle
        // φc = π / 2 - φ, which introduces a 1 ULP inconsistency between φc and
        // φ.
        if (mm != 1.0 || i != iend) {
          EXPECT_THAT(d, AlmostEquals(expected_value_d, 0, 97));
          EXPECT_THAT(j, AlmostEquals(expected_value_j, 0, 135));
        }
        ++expected_index;
      }
    }
  }
}

// Tabulated values produced with Mathematica for m close to 1 (where there were
// bugs).
TEST_F(EllipticIntegralsTest, MathematicaMNear1) {
  auto const elliptic_integrals_expected = ReadFromTabulatedData(
      SOLUTION_DIR / "numerics" / "elliptic_integrals.proto.txt");

  for (auto const& entry : elliptic_integrals_expected.entry()) {
    double const argument_n = entry.argument(0);
    double const argument_m = entry.argument(1);
    double const argument_φ = entry.argument(2);
    double const expected_value_b = entry.value(0);
    double const expected_value_d = entry.value(1);
    double const expected_value_j = entry.value(2);

    double actual_value_b;
    double actual_value_d;
    double actual_value_j;
    FukushimaEllipticBDJ(argument_φ,
                         argument_n,
                         1.0 - argument_m,
                         actual_value_b,
                         actual_value_d,
                         actual_value_j);

    EXPECT_THAT(actual_value_b, AlmostEquals(expected_value_b, 0, 7))
        << argument_n << " " << argument_m << " " << argument_φ;
    EXPECT_THAT(actual_value_d, AlmostEquals(expected_value_d, 0, 27))
        << argument_n << " " << argument_m << " " << argument_φ;
    EXPECT_THAT(actual_value_j, AlmostEquals(expected_value_j, 0, 18835))
        << argument_n << " " << argument_m << " " << argument_φ;
  }
}

}  // namespace numerics
}  // namespace principia
