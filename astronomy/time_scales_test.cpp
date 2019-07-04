
#include "astronomy/time_scales.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {
namespace internal_time_scales {

using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Nano;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using ::testing::AllOf;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;

constexpr Instant j2000_week = "1999-W52-6T12:00:00"_TT;

constexpr Instant j2000_from_tt = "2000-01-01T12:00:00"_TT;
constexpr Instant j2000_from_tai = "2000-01-01T11:59:27,816"_TAI;
constexpr Instant j2000_from_utc = "2000-01-01T11:58:55,816"_UTC;
constexpr Instant j2000_tai = "2000-01-01T12:00:00"_TAI;
constexpr Instant j2000_tai_from_tt = "2000-01-01T12:00:32,184"_TT;

class TimeScalesTest : public testing::Test {
 protected:
  static constexpr bool OneMicrosecondApart(Instant const& t1,
                                            Instant const& t2) {
    return t1 - t2 <= 1 * Micro(Second) &&
           t1 - t2 >= - 1 * Micro(Second);
  }
};

using TimeScalesDeathTest = TimeScalesTest;

// The checks are giant boolean expressions which are entirely repeated in the
// error message; we try to match the relevant part.

TEST_F(TimeScalesDeathTest, LeaplessScales) {
  EXPECT_DEATH("2015-06-30T23:59:60"_TT, "!tt.time...is_leap_second..");
  EXPECT_DEATH("2015-06-30T23:59:60"_TAI, "!tai.time...is_leap_second..");
  EXPECT_DEATH("2015-06-30T23:59:60"_UT1, "!ut1.time...is_leap_second..");
}

TEST_F(TimeScalesDeathTest, BeforeRange) {
  EXPECT_DEATH("1960-12-31T23:59:59,999"_UTC, "IsValidStretchyUTC");
  EXPECT_DEATH("1830-04-10T23:59:59,999"_UT1,
               "mjd.TimeSince20000101T120000Z.ut1.. >= "
               "experimental_eop_c02.front...ut1_mjd");
}

TEST_F(TimeScalesDeathTest, WarWasBeginning) {
  EXPECT_DEATH("2101-01-01T00:00:00"_UTC, "leap_seconds.size");
  EXPECT_DEATH("2101-01-01T00:00:00"_UT1,
               "TimeSince20000101T120000Z.ut1. < eop_c04.back...ut1..");
}

TEST_F(TimeScalesDeathTest, FirstUnknownUTC) {
  EXPECT_DEATH("2020-06-30T23:59:60"_UTC, "leap_seconds.size");
  EXPECT_DEATH("2020-06-30T24:00:00"_UTC, "leap_seconds.size");
  EXPECT_DEATH("2020-07-01T00:00:00"_UTC, "leap_seconds.size");
}

TEST_F(TimeScalesDeathTest, StretchyLeaps) {
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

TEST_F(TimeScalesDeathTest, ModernLeaps) {
  EXPECT_DEATH("2015-12-31T23:59:60"_UTC, "IsValidModernUTC");
}

TEST_F(TimeScalesTest, ConstexprJ2000) {
  static_assert(j2000_week == J2000, "");
  static_assert(j2000_from_tt == J2000, "");
  static_assert(j2000_from_tai == J2000, "");
  static_assert(j2000_from_utc == J2000, "");
  static_assert(j2000_tai == j2000_tai_from_tt, "");
  static_assert(j2000_tai - J2000 == 32.184 * Second, "");
}

TEST_F(TimeScalesTest, ConstexprWeeks) {
  // Check that week dates that go to the previous year work.
  static_assert("1914-W01-1T00:00:00"_TT == "19131229T000000"_TT, "");
}

TEST_F(TimeScalesTest, ConstexprMJD2000) {
  constexpr Instant mjd51544_utc = "2000-01-01T00:00:00"_UTC;
  constexpr Instant mjd51544_utc_from_ut1 =
      "2000-01-01T00:00:00,355"_UT1 + 388.0 * Micro(Second);

  static_assert(mjd51544_utc - mjd51544_utc_from_ut1 < 1 * Nano(Second), "");
  static_assert(mjd51544_utc - mjd51544_utc_from_ut1 > -1 * Nano(Second), "");
}

TEST_F(TimeScalesTest, ReferenceDates) {
  EXPECT_THAT("1858-11-17T00:00:00"_TT, Eq("MJD0"_TT));
  EXPECT_THAT(j2000_week, Eq(J2000));
  EXPECT_THAT(j2000_from_tt, Eq(J2000));
  EXPECT_THAT(j2000_from_tai, Eq(J2000));
  EXPECT_THAT(j2000_from_utc, Eq(J2000));
  EXPECT_THAT(j2000_tai, Eq(j2000_tai_from_tt));
  EXPECT_THAT(j2000_tai - J2000, Eq(32.184 * Second));

  // Besselian epochs.
  constexpr Instant B1900 = "1899-12-31T00:00:00"_TT + 0.8135 * Day;
  Instant const JD2415020_3135 = "JD2415020.3135"_TT;
  EXPECT_THAT(B1900, AlmostEquals(JD2415020_3135, 1));
  EXPECT_THAT(testing_utilities::AbsoluteError(JD2415020_3135, B1900),
              IsNear(0.5 * Micro(Second)));

  constexpr Instant B1950 = "1949-12-31T00:00:00"_TT + 0.9235 * Day;
  Instant const JD2433282_4235 = "JD2433282.4235"_TT;
  EXPECT_THAT(B1950, AlmostEquals(JD2433282_4235, 0));
}

TEST_F(TimeScalesTest, LeapSecond) {
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
TEST_F(TimeScalesTest, StretchyLeaps) {
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

TEST_F(TimeScalesTest, StretchyRates) {
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

TEST_F(TimeScalesTest, UT1Continuity) {
  // Continuity with TAI.  We have a fairly low resolution for UT1 at that time,
  // as well as high errors (~20 ms), and TAI was synchronized with UT2 anyway,
  // so it's not going to get much better than 100 ms.
  EXPECT_THAT(
      AbsoluteError("1958-01-01T00:00:00"_UT1, "1958-01-01T00:00:00"_TAI),
      Lt(100 * Milli(Second)));

  // Continuity at the beginning of the EOP C02 series.
  EXPECT_THAT(AbsoluteError("1961-12-31T23:59:59,000"_UT1,
                            "1961-12-31T23:59:58,967"_UTC),
              Lt(0.5 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1962-01-01T00:00:00,000"_UT1,
                            "1961-12-31T23:59:59,967"_UTC),
              Lt(0.5 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1962-01-01T00:00:00,033"_UT1,
                            "1962-01-01T00:00:00,000"_UTC),
              Lt(0.5 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1962-01-01T00:00:01,033"_UT1,
                            "1962-01-01T00:00:01,000"_UTC),
              Lt(0.5 * Milli(Second)));

  // Continuity across a stretchy UTC leap.
  EXPECT_THAT(AbsoluteError("1964-03-31T23:59:59,000"_UT1,
                            "1964-03-31T23:59:59,160"_UTC),
              Lt(0.5 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1964-03-31T23:59:59,900"_UT1,
                            "1964-03-31T23:59:60,060"_UTC),
              Lt(0.5 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1964-03-31T23:59:59,940"_UT1,
                            "1964-04-01T00:00:00,000"_UTC),
              Lt(0.5 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1964-04-01T00:00:00,000"_UT1,
                            "1964-04-01T00:00:00,060"_UTC),
              Lt(0.5 * Milli(Second)));
}

// Check the times of the lunar eclipses in LunarEclipseTest.
TEST_F(TimeScalesTest, LunarEclipses) {
  EXPECT_THAT(AbsoluteError("1950-04-02T20:44:34.0"_TT,
                            "1950-04-02T20:44:04.8"_UT1),
              Lt(14 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1950-04-02T20:49:16.7"_TT,
                            "1950-04-02T20:48:47.5"_UT1),
              Lt(14 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1950-09-26T04:17:11.4"_TT,
                            "1950-09-26T04:16:42.1"_UT1),
              Lt(86 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1950-09-26T04:21:55.5"_TT,
                            "1950-09-26T04:21:26.1"_UT1),
              Lt(15 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1951-03-23T10:37:33.2"_TT,
                            "1951-03-23T10:37:03.7"_UT1),
              Lt(92 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1951-03-23T10:50:16.8"_TT,
                            "1951-03-23T10:49:47.3"_UT1),
              Lt(92 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1951-09-15T12:27:06.3"_TT,
                            "1951-09-15T12:26:36.6"_UT1),
              Lt(99 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1951-09-15T12:38:51.5"_TT,
                            "1951-09-15T12:38:21.8"_UT1),
              Lt(99 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1952-02-11T00:28:39.9"_TT,
                            "1952-02-11T00:28:10.0"_UT1),
              Lt(69 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1952-02-11T00:39:47.6"_TT,
                            "1952-02-11T00:39:17.7"_UT1),
              Lt(69 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("1952-08-05T19:40:29.4"_TT,
                            "1952-08-05T19:39:59.3"_UT1),
              Lt(57 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("1952-08-05T19:47:54.8"_TT,
                            "1952-08-05T19:47:24.7"_UT1),
              Lt(57 * Milli(Second)));

  EXPECT_THAT(AbsoluteError("2000-01-21T04:41:30.5"_TT,
                            "2000-01-21T04:40:26.7"_UT1),
              Lt(45 * Milli(Second)));
  EXPECT_THAT(AbsoluteError("2000-01-21T04:44:34.5"_TT,
                            "2000-01-21T04:43:30.6"_UT1),
              Lt(56 * Milli(Second)));

  EXPECT_THAT("2048-01-01T06:53:54.8"_TT - "2048-01-01T06:52:23.6"_TT,
              AlmostEquals(91.2 * Second, 3e6, 4e6));
  EXPECT_THAT("2048-01-01T06:58:19.8"_TT - "2048-01-01T06:56:48.6"_TT,
              AlmostEquals(91.2 * Second, 3e6, 4e6));
}

TEST_F(TimeScalesTest, JulianDate) {
  static_assert("2010-01-04T00:00:00.108"_TT ==
                "JD2455200.50000125"_TT, "Dates differ");
  static_assert("2010-01-04T00:00:00.000"_TT ==
                "JD2455200.50000"_TT, "Dates differ");
  static_assert("2010-01-04T12:00:00.000"_TT ==
                "JD2455201.00000"_TT, "Dates differ");
  static_assert("2010-01-04T18:00:00.000"_TT ==
                "JD2455201.25000"_TT, "Dates differ");
  static_assert(OneMicrosecondApart("2010-01-04T02:57:46.659"_TT,
                                    "JD2455200.623456701388"_TT),
                "Dates differ");
  static_assert("2000-01-01T00:00:00"_TT ==
                "JD2451544.5"_TT, "Dates differ");
  static_assert("2000-01-01T12:00:00"_TT ==
                "JD2451545"_TT, "Dates differ");

  EXPECT_THAT("JD2451545"_TT, Eq(j2000_week));
  EXPECT_THAT("JD2455201.00000"_TT, Eq("2010-01-04T12:00:00.000"_TT));

  double const jd = 2457662.55467;
  Instant const date = Instant() + (jd - 2451545.0) * Day;
  EXPECT_THAT("JD2457662.55467"_TT,
              AllOf(Lt(date + 1.0 * Second), Gt(date - 1.0 * Second)))
      << date - "JD2457662.55467"_TT;
}

TEST_F(TimeScalesTest, ModifiedJulianDate) {
  static_assert("2010-01-04T00:00:00.123"_TT == "MJD55200.0000014236111"_TT,
                "Dates differ");
  static_assert("2010-01-04T00:00:00.000"_TT == "MJD55200.00000"_TT,
                "Dates differ");
  static_assert("2010-01-04T12:00:00.000"_TT == "MJD55200.50000"_TT,
                "Dates differ");
  static_assert("2010-01-04T18:00:00.000"_TT == "MJD55200.75000"_TT,
                "Dates differ");
  static_assert(OneMicrosecondApart("2010-01-04T02:57:46.659"_TT,
                                    "MJD55200.123456701388"_TT),
                "Dates differ");

  EXPECT_THAT("MJD0.0"_TT, Eq("1858-11-17T00:00:00"_TT));
  EXPECT_THAT("MJD55200.0000"_TT, Eq("2010-01-04T00:00:00.000"_TT));
  EXPECT_THAT("MJD55200.0000014236111"_TT, Eq("2010-01-04T00:00:00.123"_TT));
}

TEST_F(TimeScalesDeathTest, JulianDateUTC) {
  EXPECT_DEATH("JD2451545"_UTC, "size > 0");
  EXPECT_DEATH("MJD55200.123"_UTC, "size > 0");
}

TEST_F(TimeScalesTest, EarthRotationAngle) {
  constexpr double revolutions_at_j2000_ut1 = 0.7790572732640;
  constexpr double excess_revolutions_per_ut1_day = 0.00273781191135448;

  // Round-trip from UT1, comparing with the direct computation from UT1.
  static_assert(EarthRotationAngle("JD2451545.0"_UT1) ==
                    2 * π * Radian * revolutions_at_j2000_ut1,
                "Angles differ");
  EXPECT_THAT(
      EarthRotationAngle("JD2455200.0"_UT1),
      AlmostEquals(2 * π * Radian *
                       (revolutions_at_j2000_ut1 +
                        excess_revolutions_per_ut1_day * (2455200 - 2451545)),
                   142));
  EXPECT_THAT(
      EarthRotationAngle("JD2455200.623456701388"_UT1),
      AlmostEquals(
          2 * π * Radian *
              (0.623456701388 + revolutions_at_j2000_ut1 +
               excess_revolutions_per_ut1_day * (5200.623456701388 - 1545) - 1),
          134));

  // Compare with the WGCCRE 2009 elements.
  EXPECT_THAT(
      (EarthRotationAngle(J2000) - π / 2 * Radian) - 190.147 * Degree,
      IsNear(0.0469 * Degree));
  EXPECT_THAT((EarthRotationAngle("2000-01-01T23:00:00"_TT) -
               EarthRotationAngle("2000-01-01T01:00:00"_TT)) /
                      (22 * Hour) -
                  360.9856235 * (Degree / Day),
              IsNear(-0.0000149 * Degree / Day));
  EXPECT_THAT((EarthRotationAngle("2010-01-01T23:00:00"_TT) -
               EarthRotationAngle("2010-01-01T01:00:00"_TT)) /
                      (22 * Hour) -
                  360.9856235 * (Degree / Day),
              IsNear(-0.0000137 * Degree / Day));
}

TEST_F(TimeScalesTest, GNSS) {
  // BeiDou Navigation Satellite System
  // Signal In Space Interface Control Document
  // Open Service Signals B1C and B2a (Test Version),
  // 3.3 Time System.
  // The start epoch of BDT is 00:00:00 on January 1, 2006 of Coordinated
  // Universal Time (UTC).
  EXPECT_THAT("2006-01-01T00:00:00"_北斗, Eq("2006-01-01T00:00:00"_UTC));

  // Galileo OS SIS ICD, Issue 1.1.
  // 5.1.2. Galileo System Time (GST).
  // The GST start epoch shall be 00:00 UT on Sunday 22nd August 1999 (midnight
  // between 21st and 22nd August). At the start epoch, GST shall be ahead of
  // UTC by thirteen (13) leap seconds.
  EXPECT_THAT("1999-08-22T00:00:13"_GPS, Eq("1999-08-22T00:00:00"_UTC));
}

}  // namespace internal_time_scales
}  // namespace astronomy
}  // namespace principia
