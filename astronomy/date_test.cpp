
#include "astronomy/date.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {
namespace internal_date {

using quantities::si::Micro;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Lt;

class DateTest : public testing::Test {};

using DateDeathTest = DateTest;

// The checks are giant boolean expressions which are entirely repeated in the
// error message; we try to match the relevant part.

#if !((PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL) && WE_LIKE_N3599)

TEST_F(DateDeathTest, InvalidCalendarDate) {
  EXPECT_DEATH("2001-04-00T12:00:00Z"_TT, "day >= 1");
  EXPECT_DEATH("2001-02-29T12:00:00Z"_TT, "day <= month_length");
  EXPECT_DEATH("2001-03-32T12:00:00Z"_TT, "day <= month_length");
  EXPECT_DEATH("2001-04-31T12:00:00Z"_TT, "day <= month_length");
  EXPECT_DEATH("2001-00-01T12:00:00Z"_TT, "month >= 1");
  EXPECT_DEATH("2001-13-01T12:00:00Z"_TT, "month <= 12");
  EXPECT_DEATH("2001-00-01T12:00:00Z"_TT, "month >= 1");
  EXPECT_DEATH("1582-01-01T12:00:00Z"_TT, "year >= 1583");
}

TEST_F(DateDeathTest, InvalidTime) {
  EXPECT_DEATH("2001-01-01T25:00:00Z"_TT, "hour_ <= 23");
  EXPECT_DEATH("2001-01-01T24:01:00Z"_TT, "minute_ == 0");
  EXPECT_DEATH("2001-01-01T24:00:01Z"_TT, "second_ == 0");
  EXPECT_DEATH("2001-01-01T00:60:00Z"_TT, "minute_ <= 59");
  EXPECT_DEATH("2001-01-01T00:00:60Z"_TT, "second_ <= 59");
  EXPECT_DEATH("2001-01-01T23:59:61Z"_TT, "second_ == 60");
}

TEST_F(DateDeathTest, InvalidDateTime) {
  EXPECT_DEATH("2001-01-01T23:59:60Z"_TT, "");
}

#endif

namespace {

constexpr Instant j2000_week = "1999-W52-6T12:00:00Z"_TT;

constexpr Instant j2000_from_tt = "2000-01-01T12:00:00Z"_TT;
constexpr Instant j2000_from_tai = "2000-01-01T11:59:27,816Z"_TAI;
constexpr Instant j2000_from_utc = "2000-01-01T11:58:55,816Z"_UTC;
constexpr Instant j2000_tai = "2000-01-01T12:00:00Z"_TAI;
constexpr Instant j2000_tai_from_tt = "2000-01-01T12:00:32,184Z"_TT;

static_assert(j2000_week == J2000, "");
static_assert(j2000_from_tt == J2000, "");
static_assert(j2000_from_tai == J2000, "");
static_assert(j2000_from_utc == J2000, "");
static_assert(j2000_tai == j2000_tai_from_tt, "");
static_assert(j2000_tai - J2000 == 32.184 * Second, "");

// Check that week dates that go to the previous year work.
static_assert("1914-W01-1T00:00:00Z"_TT == "19131229T000000Z"_TT, "");

}  // namespace

TEST_F(DateTest, ReferenceDates) {
  EXPECT_THAT("1858-11-17T00:00:00Z"_TT, Eq(ModifiedJulianDate(0)));
  EXPECT_THAT(j2000_week, Eq(J2000));
  EXPECT_THAT(j2000_from_tt, Eq(J2000));
  EXPECT_THAT(j2000_from_tai, Eq(J2000));
  EXPECT_THAT(j2000_from_utc, Eq(J2000));
  EXPECT_THAT(j2000_tai, Eq(j2000_tai_from_tt));
  EXPECT_THAT(j2000_tai - J2000, Eq(32.184 * Second));

  // Besselian epochs.
  constexpr Instant B1900 = "1899-12-31T00:00:00Z"_TT + 0.8135 * Day;
  Instant const JD2415020_3135 = JulianDate(2415020.3135);
  EXPECT_THAT(B1900, AlmostEquals(JD2415020_3135, 51));
  EXPECT_THAT(testing_utilities::AbsoluteError(JD2415020_3135, B1900),
              AllOf(Ge(10 * Micro(Second)), Lt(100 * Micro(Second))));

  constexpr Instant B1950 = "1949-12-31T00:00:00Z"_TT + 0.9235 * Day;
  Instant const JD2433282_4235 = JulianDate(2433282.4235);
  EXPECT_THAT(B1950, AlmostEquals(JD2433282_4235, 26));
  EXPECT_THAT(testing_utilities::AbsoluteError(JD2433282_4235, B1950),
              AllOf(Ge(1 * Micro(Second)), Lt(10 * Micro(Second))));
}

TEST_F(DateTest, LeapSecond) {
  Instant const eleven_fifty_nine_and_fifty_eight_seconds =
      "2015-06-30T23:59:58Z"_UTC;
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 1 * Second,
              Eq("2015-06-30T23:59:59Z"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 1.5 * Second,
              Eq("2015-06-30T23:59:59,5Z"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 2 * Second,
              Eq("2015-06-30T23:59:60Z"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 2.5 * Second,
              Eq("2015-06-30T23:59:60,5Z"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 3 * Second,
              Eq("2015-06-30T24:00:00Z"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 3 * Second,
              Eq("2015-07-01T00:00:00Z"_UTC));
}

}  // namespace internal_date
}  // namespace astronomy
}  // namespace principia
