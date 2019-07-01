
#include "numerics/elliptic_functions.hpp"

#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/elliptic_integrals.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace numerics {

using quantities::Angle;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::ReadFromTabulatedData;
using ::testing::Le;

class EllipticFunctionsTest : public ::testing::Test {};

// The test values found in Fukushima's xgscd.txt file.
TEST_F(EllipticFunctionsTest, Xgscd) {
  auto const xgscd_expected =
      ReadFromTabulatedData(SOLUTION_DIR / "numerics" / "xgscd.proto.txt");
  double Δmc, mc, m, du, u, s, c, d;
  int jend, iend;
  jend = 10;
  iend = 8;
  Δmc = 1.0 / static_cast<double>(jend);
  std::printf("%10s%10s%25s%25s%25s\n", "m", "u", "s", "c", "d");
  int expected_index = 0;
  for (int j = 1; j <= jend; ++j) {
    mc = static_cast<double>(j) * Δmc;
    m = 1.0 - mc;
    du = EllipticK(mc) / static_cast<double>(iend);
    for (int i = 0; i <= iend * 8; ++i) {
      u = du * static_cast<double>(i);
      JacobiSNCNDN(u, mc, s, c, d);
      std::printf("%10.5f%10.5f%25.15e%25.15e%25.15e\n", m, u, s, c, d);

      auto const& expected_entry = xgscd_expected.entry(expected_index);
      auto const expected_argument_m = expected_entry.argument(0);
      auto const expected_argument_u = expected_entry.argument(1);
      auto const expected_value_s = expected_entry.value(0);
      auto const expected_value_c = expected_entry.value(1);
      auto const expected_value_d = expected_entry.value(2);
      EXPECT_THAT(m, IsNear(expected_argument_m, 1.001));
      EXPECT_THAT(u, IsNear(expected_argument_u, 1.001));
      EXPECT_THAT(s, AlmostEquals(expected_value_s, 0, 2));
      EXPECT_THAT(c, AlmostEquals(expected_value_c, 0, 3));
      EXPECT_THAT(d, AlmostEquals(expected_value_d, 0, 1));
      ++expected_index;
    }
    std::printf("\n");
  }
}

TEST_F(EllipticFunctionsTest, Mathematica) {
  auto const elliptic_functions_expected = ReadFromTabulatedData(
      SOLUTION_DIR / "numerics" / "elliptic_functions.proto.txt");

  for (auto const& entry : elliptic_functions_expected.entry()) {
    double const argument_u = entry.argument(0);
    double const argument_m = entry.argument(1);
    double const expected_value_s = entry.value(0);
    double const expected_value_c = entry.value(1);
    double const expected_value_d = entry.value(2);
    Angle const expected_value_a = entry.value(3) * Radian;

    double actual_value_s;
    double actual_value_c;
    double actual_value_d;
    JacobiSNCNDN(argument_u,
                 1.0 - argument_m,
                 actual_value_s,
                 actual_value_c,
                 actual_value_d);

    Angle const actual_value_a = JacobiAmplitude(argument_u, 1.0 - argument_m);

    // TODO(phl): The errors are uncomfortably large here.  Figure out what's
    // happening.
    EXPECT_THAT(actual_value_s, AlmostEquals(expected_value_s, 0, 12507))
        << argument_u << " " << argument_m;
    EXPECT_THAT(actual_value_c, AlmostEquals(expected_value_c, 0, 7648))
        << argument_u << " " << argument_m;
    EXPECT_THAT(actual_value_d, AlmostEquals(expected_value_d, 0, 85))
        << argument_u << " " << argument_m;

    EXPECT_THAT(actual_value_a, AlmostEquals(expected_value_a, 0, 22073))
        << argument_u << " " << argument_m;
  }
}

TEST_F(EllipticFunctionsTest, Monotonicity) {
  for (double const mc : {0.01, 0.1, 0.5}) {
    double const k = EllipticK(mc);
    for (int i = -5; i <= 5; ++i) {
      double const u = 2.0 * i * k;
      Angle const am = JacobiAmplitude(u, mc);
      LOG(ERROR) << "Testing around " << u << " and mc = " << mc
                 << " starting with am = " << am << " (i = " << i << ")";
      double u₊ = u;
      double u₋ = u;
      Angle previous_am₊ = am;
      Angle previous_am₋ = am;
      for (int j = 0; j < 1'000; ++j) {
        u₊ = std::nextafter(u₊, std::numeric_limits<double>::infinity());
        u₋ = std::nextafter(u₋, -std::numeric_limits<double>::infinity());
        Angle const am₊ = JacobiAmplitude(u₊, mc);
        Angle const am₋ = JacobiAmplitude(u₋, mc);
        ASSERT_THAT(previous_am₊, Le(am₊)) << u₊ << " " << mc << " " << j;
        ASSERT_THAT(am₋, Le(previous_am₋)) << u₋ << " " << mc << " " << j;
        previous_am₋ = am₋;
        previous_am₊ = am₊;
      }
    }
  }
}

}  // namespace numerics
}  // namespace principia
