
#include "astronomy/date.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {
namespace internal_date {

using quantities::si::Micro;
using quantities::si::Milli;
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

#if !((PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL) && \
      WE_LIKE_N3599)

TEST_F(DateDeathTest, InvalidCalendarDate) {
  EXPECT_DEATH("2001-04-00T12:00:00"_TT, "day >= 1");
  EXPECT_DEATH("2001-02-29T12:00:00"_TT, "day <= month_length");
  EXPECT_DEATH("2001-03-32T12:00:00"_TT, "day <= month_length");
  EXPECT_DEATH("2001-04-31T12:00:00"_TT, "day <= month_length");
  EXPECT_DEATH("2001-00-01T12:00:00"_TT, "month >= 1");
  EXPECT_DEATH("2001-13-01T12:00:00"_TT, "month <= 12");
  EXPECT_DEATH("2001-00-01T12:00:00"_TT, "month >= 1");
  EXPECT_DEATH("1582-01-01T12:00:00"_TT, "year >= 1583");
}

TEST_F(DateDeathTest, InvalidTime) {
  EXPECT_DEATH("2001-01-01T25:00:00"_TT, "hour_ <= 23");
  EXPECT_DEATH("2001-01-01T24:01:00"_TT, "minute_ == 0");
  EXPECT_DEATH("2001-01-01T24:00:01"_TT, "second_ == 0");
  EXPECT_DEATH("2001-01-01T00:60:00"_TT, "minute_ <= 59");
  EXPECT_DEATH("2001-01-01T00:00:60"_TT, "second_ <= 59");
  EXPECT_DEATH("2001-01-01T23:59:61"_TT, "second_ == 60");
}

TEST_F(DateDeathTest, InvalidDateTime) {
  EXPECT_DEATH("2001-01-01T23:59:60"_TT, "");
}

TEST_F(DateDeathTest, StretchyLeaps) {
  EXPECT_DEATH("1961-07-31T23:59:59,950"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1963-10-31T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1964-03-31T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1964-08-31T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1964-12-31T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1965-02-28T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1965-06-30T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1965-08-31T23:59:60,100"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1968-01-31T23:59:59,900"_UTC, "IsValidStretchyUTC");
}

#endif

namespace {

constexpr Instant j2000_week = "1999-W52-6T12:00:00"_TT;

constexpr Instant j2000_from_tt = "2000-01-01T12:00:00"_TT;
constexpr Instant j2000_from_tai = "2000-01-01T11:59:27,816"_TAI;
constexpr Instant j2000_from_utc = "2000-01-01T11:58:55,816"_UTC;
constexpr Instant j2000_tai = "2000-01-01T12:00:00"_TAI;
constexpr Instant j2000_tai_from_tt = "2000-01-01T12:00:32,184"_TT;

static_assert(j2000_week == J2000, "");
static_assert(j2000_from_tt == J2000, "");
static_assert(j2000_from_tai == J2000, "");
static_assert(j2000_from_utc == J2000, "");
static_assert(j2000_tai == j2000_tai_from_tt, "");
static_assert(j2000_tai - J2000 == 32.184 * Second, "");

// Check that week dates that go to the previous year work.
static_assert("1914-W01-1T00:00:00"_TT == "19131229T000000"_TT, "");

}  // namespace

TEST_F(DateTest, ReferenceDates) {
  EXPECT_THAT("1858-11-17T00:00:00"_TT, Eq(ModifiedJulianDate(0)));
  EXPECT_THAT(j2000_week, Eq(J2000));
  EXPECT_THAT(j2000_from_tt, Eq(J2000));
  EXPECT_THAT(j2000_from_tai, Eq(J2000));
  EXPECT_THAT(j2000_from_utc, Eq(J2000));
  EXPECT_THAT(j2000_tai, Eq(j2000_tai_from_tt));
  EXPECT_THAT(j2000_tai - J2000, Eq(32.184 * Second));

  // Besselian epochs.
  constexpr Instant B1900 = "1899-12-31T00:00:00"_TT + 0.8135 * Day;
  Instant const JD2415020_3135 = JulianDate(2415020.3135);
  EXPECT_THAT(B1900, AlmostEquals(JD2415020_3135, 51));
  EXPECT_THAT(testing_utilities::AbsoluteError(JD2415020_3135, B1900),
              AllOf(Ge(10 * Micro(Second)), Lt(100 * Micro(Second))));

  constexpr Instant B1950 = "1949-12-31T00:00:00"_TT + 0.9235 * Day;
  Instant const JD2433282_4235 = JulianDate(2433282.4235);
  EXPECT_THAT(B1950, AlmostEquals(JD2433282_4235, 26));
  EXPECT_THAT(testing_utilities::AbsoluteError(JD2433282_4235, B1950),
              AllOf(Ge(1 * Micro(Second)), Lt(10 * Micro(Second))));
}

TEST_F(DateTest, LeapSecond) {
  Instant const eleven_fifty_nine_and_fifty_eight_seconds =
      "2015-06-30T23:59:58"_UTC;
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 1 * Second,
              Eq("2015-06-30T23:59:59"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 1.5 * Second,
              Eq("2015-06-30T23:59:59,5"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 2 * Second,
              Eq("2015-06-30T23:59:60"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 2.5 * Second,
              Eq("2015-06-30T23:59:60,5"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 3 * Second,
              Eq("2015-06-30T24:00:00"_UTC));
  EXPECT_THAT(eleven_fifty_nine_and_fifty_eight_seconds + 3 * Second,
              Eq("2015-07-01T00:00:00"_UTC));

  constexpr Instant end_of_december_2016 = "2016-12-31T24:00:00"_UTC;
  EXPECT_THAT(end_of_december_2016 - 2 * Second, Eq("2016-12-31T23:59:59"_UTC));
  EXPECT_THAT(end_of_december_2016 - 1 * Second, Eq("2016-12-31T23:59:60"_UTC));
  EXPECT_THAT(end_of_december_2016 - 0 * Second, Eq("2017-01-01T00:00:00"_UTC));
  EXPECT_THAT("2016-12-31T23:59:59"_UTC - "2016-12-31T23:59:59"_TAI,
              Eq(36 * Second));
  EXPECT_THAT("2017-01-01T00:00:00"_UTC - "2017-01-01T00:00:00"_TAI,
              Eq(37 * Second));
}

// See the list of steps at
// https://hpiers.obspm.fr/iers/bul/bulc/TimeSteps.history.
// Note that while the same file is used to check that the date string is valid
// with respect to positive or negative leap seconds, the actual conversion is
// based exclusively on https://hpiers.obspm.fr/iers/bul/bulc/UTC-TAI.history,
// so this provides some sort of cross-checking.
TEST_F(DateTest, StretchyLeaps) {
  EXPECT_THAT(AbsoluteError("1961-07-31T24:00:00,000"_UTC - 0.050 * Second,
                            "1961-07-31T23:59:59,900"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1961-08-01T00:00:00"_UTC, "1961-08-01T00:00:01,648"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1963-10-31T24:00:00,000"_UTC - 0.100 * Second,
                            "1963-10-31T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1963-11-01T00:00:00"_UTC, "1963-11-01T00:00:02,697"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1964-03-31T24:00:00,000"_UTC - 0.100 * Second,
                            "1964-03-31T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1964-04-01T00:00:00"_UTC, "1964-04-01T00:00:02,984"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1964-08-31T24:00:00,000"_UTC - 0.100 * Second,
                            "1964-08-31T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1964-09-01T00:00:00"_UTC, "1964-09-01T00:00:03,282"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1964-12-31T24:00:00,000"_UTC - 0.100 * Second,
                            "1964-12-31T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1965-01-01T00:00:00"_UTC, "1965-01-01T00:00:03,540"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1965-02-28T24:00:00,000"_UTC - 0.100 * Second,
                            "1965-02-28T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1965-03-01T00:00:00"_UTC, "1965-03-01T00:00:03,717"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1965-06-30T24:00:00,000"_UTC - 0.100 * Second,
                            "1965-06-30T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1965-07-01T00:00:00"_UTC, "1965-07-01T00:00:03,975"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1965-08-31T24:00:00,000"_UTC - 0.100 * Second,
                            "1965-08-31T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1965-09-01T00:00:00"_UTC, "1965-09-01T00:00:04,155"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1968-01-31T24:00:00,000"_UTC - 0.100 * Second,
                            "1968-01-31T23:59:59,800"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1968-02-01T00:00:00"_UTC, "1968-02-01T00:00:06,186"_TAI),
      Lt(0.5 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1971-12-31T24:00:00,000"_UTC - 0.107'7580 * Second,
                            "1971-12-31T23:59:60,000"_UTC),
              Lt(1 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError("1972-01-01T00:00:00"_UTC, "1972-01-01T00:00:10,000"_TAI),
      Lt(0.5 * Milli(Second)));
}

TEST_F(DateTest, StretchyRates) {
  // Check that cancellations aren't destroying the test.
  EXPECT_NE("1961-01-01T00:00:00"_UTC + 1 * Minute / (1 - 150e-10),
            "1961-01-01T00:00:00"_UTC + 1 * Minute / (1 - 130e-10));

  quantities::Time utc_minute;
  utc_minute = 1 * Minute / (1 - 150e-10);

  EXPECT_THAT("1961-01-01T00:00:00"_UTC + utc_minute,
              Eq("1961-01-01T00:01:00"_UTC));
  EXPECT_THAT("1961-12-31T23:59:00"_UTC + utc_minute,
              Eq("1961-12-31T24:00:00"_UTC));

  utc_minute = 1 * Minute / (1 - 130e-10);
  EXPECT_THAT("1962-01-01T00:00:00"_UTC + utc_minute,
              Eq("1962-01-01T00:01:00"_UTC));
  EXPECT_THAT("1963-12-31T23:59:00"_UTC + utc_minute,
              Eq("1963-12-31T24:00:00"_UTC));

  utc_minute = 1 * Minute / (1 - 150e-10);
  EXPECT_THAT("1964-01-01T00:00:00"_UTC + utc_minute,
              Eq("1964-01-01T00:01:00"_UTC));
  EXPECT_THAT("1965-12-31T23:59:00"_UTC + utc_minute,
              AlmostEquals("1965-12-31T24:00:00"_UTC, 1));

  utc_minute = 1 * Minute / (1 - 300e-10);
  EXPECT_THAT("1966-01-01T00:00:00"_UTC + utc_minute,
              Eq("1966-01-01T00:01:00"_UTC));
  EXPECT_THAT("1971-12-31T23:58:00"_UTC + utc_minute,
              Eq("1971-12-31T23:59:00"_UTC));

  utc_minute = 1 * Minute;
  EXPECT_THAT("1972-01-01T00:00:00"_UTC + utc_minute,
              Eq("1972-01-01T00:01:00"_UTC));
  EXPECT_THAT("2000-01-01T00:00:00"_UTC + utc_minute,
              Eq("2000-01-01T00:01:00"_UTC));
}

}  // namespace internal_date
}  // namespace astronomy
}  // namespace principia
