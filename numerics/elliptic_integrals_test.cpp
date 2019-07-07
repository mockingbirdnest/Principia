
#include <filesystem>
#include <vector>
#include <utility>

#include "numerics/elliptic_integrals.hpp"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Angle;
using quantities::si::Radian;
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

  constexpr int lend = 5;
  constexpr int kend = 5;
  constexpr int iend = 4;
  constexpr double Δnc = 1.0 / (lend - 1);
  constexpr double Δmc = 1.0 / (kend - 1);
  constexpr Angle Δφ = (π / 2) * Radian / iend;
  std::printf(
      "%10s%10s%10s%25s%25s%25s\n", "n", "m", "φ / π", "elb", "eld", "elj");

  int expected_index = 0;
  for (int l = 1; l <= lend; ++l) {
    double nc = (l - 1) * Δnc;
    if (nc <= 2.44e-4) {
      nc = 2.44e-4;
    }
    double const nn = 1.0 - nc;
    for (int k = 1; k <= kend; ++k) {
      std::printf("\n");
      double mc = (k - 1) * Δmc;
      if (mc <= 0.0) {
        mc = 1.21e-32;
      }
      double const mm = 1.0 - mc;
      for (int i = 0; i <= iend; ++i) {
        double b, d, j;
        Angle const φ = Δφ * i;

        // The following is useful for tracing the actual machine numbers passed
        // to the elliptic integrals close to the pole.  We want to ensure that
        // the parameters have the exact value used in Mathematica computations.
        if (mm == 1.0 && i == iend) {
          EXPECT_THAT(mc, AlmostEquals(1.21000'00000'00000'04719e-32, 0));
          EXPECT_THAT(φ, AlmostEquals(1.57079'63267'94896'55800 * Radian, 0));
          LOG(INFO) << " mc = " << quantities::DebugString(mc, 30) << u8" φ = "
                    << quantities::DebugString(φ, 30);
          if (l == 1) {
            EXPECT_THAT(nn, AlmostEquals(9.99755'99999'99999'78023e-1, 0));
            LOG(INFO) << " n = " << quantities::DebugString(nn, 30);
          }
        }

        FukushimaEllipticBDJ(φ, nn, mc, b, d, j);
        std::printf("%10.5f%10.5f%10.5f%25.15f%25.15f%25.15f\n",
                    nn,
                    mm,
                    φ / (π * Radian),
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
        EXPECT_THAT(φ / (π * Radian),
                    IsNear(expected_argument_φ_over_π, 1.001));

        // The relatively large errors on D and J below are not significant.
        // They come from the fact that the data given by Fukushima has 15
        // digits after the decimal point (and no exponent) but D and J can get
        // quite small on this dataset (around 0.2) so the given expected values
        // don't have 53 bits in their mantissa.
        EXPECT_THAT(b, AlmostEquals(expected_value_b, 0, 8));
        EXPECT_THAT(d, AlmostEquals(expected_value_d, 0, 97));
        EXPECT_THAT(j, AlmostEquals(expected_value_j, 0, 135));
        ++expected_index;
      }
    }
  }
}

TEST_F(EllipticIntegralsTest, MathematicaBivariate) {
  auto const bivariate_elliptic_integrals_expected = ReadFromTabulatedData(
      SOLUTION_DIR / "numerics" / "bivariate_elliptic_integrals.proto.txt");

  for (auto const& entry : bivariate_elliptic_integrals_expected.entry()) {
    Angle const argument_φ = entry.argument(0) * Radian;
    double const argument_m = entry.argument(1);
    double const expected_value_b = entry.value(0);
    double const expected_value_d = entry.value(1);
    double const expected_value_e = entry.value(2);
    double const expected_value_f = entry.value(3);

    double actual_value_b;
    double actual_value_d;
    double actual_value_j;  // Ignored.
    FukushimaEllipticBDJ(argument_φ,
                         /*n=*/1,
                         1.0 - argument_m,
                         actual_value_b,
                         actual_value_d,
                         actual_value_j);

    double actual_value_e;
    double actual_value_f;
    double actual_value_ᴨ;  // Ignored.
    EllipticEFΠ(argument_φ,
                /*n=*/1,
                1.0 - argument_m,
                actual_value_e,
                actual_value_f,
                actual_value_ᴨ);

    EXPECT_THAT(actual_value_b, AlmostEquals(expected_value_b, 0, 6))
        << argument_φ << " " << argument_m;
    EXPECT_THAT(actual_value_d, AlmostEquals(expected_value_d, 0, 44))
        << argument_φ << " " << argument_m;
    EXPECT_THAT(actual_value_e, AlmostEquals(expected_value_e, 0, 6))
        << argument_φ << " " << argument_m;
    EXPECT_THAT(actual_value_f, AlmostEquals(expected_value_f, 0, 6))
        << argument_φ << " " << argument_m;
  }
}

TEST_F(EllipticIntegralsTest, MathematicaTrivariate) {
  auto const trivariate_elliptic_integrals_expected = ReadFromTabulatedData(
      SOLUTION_DIR / "numerics" / "trivariate_elliptic_integrals.proto.txt");

  for (auto const& entry : trivariate_elliptic_integrals_expected.entry()) {
    Angle const argument_φ = entry.argument(0) * Radian;
    double const argument_n = entry.argument(1);
    double const argument_m = entry.argument(2);
    double const expected_value_j = entry.value(0);
    double const expected_value_ᴨ = entry.value(1);

    double actual_value_b;  // Ignored.
    double actual_value_d;  // Ignored.
    double actual_value_j;
    FukushimaEllipticBDJ(argument_φ,
                         argument_n,
                         1.0 - argument_m,
                         actual_value_b,
                         actual_value_d,
                         actual_value_j);

    double actual_value_e;  // Ignored.
    double actual_value_f;  // Ignored.
    double actual_value_ᴨ;
    EllipticEFΠ(argument_φ,
                argument_n,
                1.0 - argument_m,
                actual_value_e,
                actual_value_f,
                actual_value_ᴨ);

    // TODO(phl): The error is uncomfortably large here.  Figure out what's
    // happening.
    EXPECT_THAT(actual_value_j, AlmostEquals(expected_value_j, 0, 25004))
        << argument_φ << " " << argument_n << " " << argument_m;
    EXPECT_THAT(actual_value_ᴨ, AlmostEquals(expected_value_ᴨ, 0, 3934))
        << argument_φ << " " << argument_n << " " << argument_m;
  }
}

}  // namespace numerics
}  // namespace principia
