#include "geometry/instant.hpp"

#include <strstream>

#include "astronomy/time_scales.hpp"
#include "geometry/sign.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace principia {
namespace geometry {

using ::testing::Eq;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_time_scales;
using namespace principia::geometry::_instant;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class InstantOutputTest : public ::testing::Test {};

TEST_F(InstantOutputTest, UniversalTime) {
  // Note that since we are roughly twelve hours from J2000, we expect
  // log₁₀(12 h / 1 s), between four and five, of the printed digits to be
  // unrepresentable: our granularity is about seven picoseconds. This is
  // obvious in the case of UTC, which is offset from TT(TAI) by an integer
  // number of milliseconds at any time.
  EXPECT_THAT((std::stringstream() << "2000-01-01T00:00:00"_UTC).str(),
              Eq("2000-01-01T00:01:04,1840000000011059 (TT)"));
  EXPECT_THAT(
      (std::stringstream() << JustAfter("2000-01-01T00:00:00"_UTC)).str(),
      Eq("2000-01-01T00:01:04,1840000000083819 (TT)"));

  EXPECT_THAT((std::stringstream() << "2000-01-01T00:00:00"_UT1).str(),
              Eq("2000-01-01T00:01:03,8286119957119809 (TT)"));

  // A second after J2000, the granularity is a couple tenths of femtoseconds.
  EXPECT_THAT((std::stringstream() << "2000-01-01T11:58:56,816"_UTC).str(),
              Eq("2000-01-01T12:00:01,0000000000000000 (TT)"));
  EXPECT_THAT(
      (std::stringstream() << JustAfter("2000-01-01T11:58:56,816"_UTC)).str(),
      Eq("2000-01-01T12:00:01,0000000000000002 (TT)"));
  // TODO(egg): This is horribly misrounded because of the 32.184 in |FromTAI|;
  // we should be doing that in |DateTime| arithmetic.
  EXPECT_THAT((std::stringstream() << "2000-01-01T11:58:56,817"_UTC).str(),
              Eq("2000-01-01T12:00:01,0009999999999977 (TT)"));

  // Times from VA501, see https://esamultimedia.esa.int/docs/esa-x-1819eng.pdf.
  // The granularity around that time is about fifteen nanoseconds.
  constexpr Instant H₀ = "1996-06-04T12:33:59"_UTC;
  constexpr Instant backup_sri_malfunction = H₀ + 36.672 * Second;
  constexpr Instant nominal_sri_malfunction = H₀ + 36.749 * Second;
  EXPECT_THAT((std::stringstream() << H₀).str(),
              Eq("1996-06-04T12:35:01,1840000003576279 (TT)"));
  EXPECT_THAT((std::stringstream() << backup_sri_malfunction).str(),
              Eq("1996-06-04T12:35:37,8560000061988831 (TT)"));
  EXPECT_THAT((std::stringstream() << nominal_sri_malfunction).str(),
              Eq("1996-06-04T12:35:37,9329999983310699 (TT)"));
}

TEST_F(InstantOutputTest, J2000) {
  EXPECT_THAT((std::stringstream() << "2000-01-01T11:59:59"_TT).str(),
              Eq("2000-01-01T11:59:59,0000000000000000 (TT)"));
  EXPECT_THAT(
      (std::stringstream() << JustAfter("2000-01-01T11:59:59"_TT)).str(),
      Eq("J2000-9.99999999999999889e-01 s (TT)"));
  EXPECT_THAT((std::stringstream() << JustBefore(J2000)).str(),
              Eq("J2000-4.94065645841246544e-324 s (TT)"));
  EXPECT_THAT((std::stringstream() << JustAfter(JustBefore(J2000))).str(),
              Eq("J2000-0.00000000000000000e+00 s (TT)"));
  EXPECT_THAT((std::stringstream() << J2000).str(),
              Eq("J2000+0.00000000000000000e+00 s (TT)"));
  EXPECT_THAT((std::stringstream() << JustAfter(J2000)).str(),
              Eq("J2000+4.94065645841246544e-324 s (TT)"));
  EXPECT_THAT(
      (std::stringstream() << JustBefore("2000-01-01T12:00:01"_TT)).str(),
      Eq("J2000+9.99999999999999889e-01 s (TT)"));
  EXPECT_THAT((std::stringstream() << "2000-01-01T12:00:01"_TT).str(),
              Eq("2000-01-01T12:00:01,0000000000000000 (TT)"));
  // 2000-01-01T11:58:55,816 UTC is J2000, and UT1 is within a second of UTC, so
  // this is within a second of J2000.
  EXPECT_THAT((std::stringstream() << "2000-01-01T11:58:55,816"_UT1).str(),
              Eq("J2000-3.54947059731458125e-01 s (TT)"));
}

TEST_F(InstantOutputTest, DistantPast) {
  // Note that near JD0, our granularity is about thirty microseconds.
  EXPECT_THAT((std::stringstream() << JustAfter("JD0.5"_TT)).str(),
              Eq("J-4712-01-02T00:00:00,0000305175781250 (TT)"));
  EXPECT_THAT((std::stringstream() << "JD0.5"_TT).str(),
              Eq("J-4712-01-02T00:00:00,0000000000000000 (TT)"));
  EXPECT_THAT((std::stringstream() << JustBefore("JD0.5"_TT)).str(),
              Eq("J2000-2.11813444800000031e+11 s (TT)"));
}

TEST_F(InstantOutputTest, DistantFuture) {
  // Likewise in the year 10000, our granularity is about thirty microseconds.
  EXPECT_THAT(
      (std::stringstream() << JustBefore("9999-12-31T24:00:00"_TT)).str(),
      Eq("9999-12-31T23:59:59,9999694824218750 (TT)"));
  EXPECT_THAT((std::stringstream() << "9999-12-31T24:00:00"_TT).str(),
              Eq("J2000+2.52455572800000000e+11 s (TT)"));
}

TEST_F(InstantOutputTest, InfinitePast) {
  EXPECT_THAT((std::stringstream() << JustAfter(InfinitePast)).str(),
              Eq("J2000-1.79769313486231571e+308 s (TT)"));
  EXPECT_THAT((std::stringstream() << InfinitePast).str(),
              Eq("J2000-inf s (TT)"));
}

TEST_F(InstantOutputTest, InfiniteFuture) {
  EXPECT_THAT((std::stringstream() << JustBefore(InfiniteFuture)).str(),
              Eq("J2000+1.79769313486231571e+308 s (TT)"));
  EXPECT_THAT((std::stringstream() << InfiniteFuture).str(),
              Eq("J2000+inf s (TT)"));
}

}  // namespace geometry
}  // namespace principia
