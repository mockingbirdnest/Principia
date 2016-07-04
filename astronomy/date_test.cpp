
#include "astronomy/date.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace astronomy {
namespace internal_date {

using ::testing::Eq;

class DateTest : public testing::Test {};

using DateDeathTest = DateTest;

// The checks are giant boolean expressions which are entirely repeated in the
// error message; we try to match the relevant part.

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

constexpr Instant j2000_week = "1999-W52-6T12:00:00Z"_TT;

constexpr Instant j2000_from_tt = "2000-01-01T12:00:00Z"_TT;
constexpr Instant j2000_from_tai = "2000-01-01T11:59:27,816Z"_TAI;
constexpr Instant j2000_from_utc = "2000-01-01T11:58:55,816Z"_UTC;
constexpr Instant j2000_tai = "2000-01-01T12:00:00Z"_TAI;
constexpr Instant j2000_tai_from_tt = "2000-01-01T12:00:32,184Z"_TT;

TEST_F(DateDeathTest, ReferenceDates) {
  EXPECT_THAT("1858-11-17T00:00:00Z"_TT, Eq(ModifiedJulianDate(0)));
  EXPECT_THAT(j2000_week, Eq(J2000));
  EXPECT_THAT(j2000_from_tt, Eq(J2000));
  EXPECT_THAT(j2000_from_tai, Eq(J2000));
  EXPECT_THAT(j2000_from_utc, Eq(J2000));
  EXPECT_THAT(j2000_tai, Eq(j2000_tai_from_tt));
  EXPECT_THAT(j2000_tai - J2000, Eq(32.184 * Second));

  // Besselian epochs.
  constexpr Instant B1900 = "1899-12-31T00:00:00Z"_TT + 0.8135 * Day;
  EXPECT_THAT(B1900, Eq(JulianDate(2415020.3135)));

  constexpr Instant B1950 = "1949-12-31T00:00:00Z"_TT + 0.9235 * Day;
  EXPECT_THAT(B1950, Eq(JulianDate(2433282.4235)));
}

}  // namespace internal_date
}  // namespace astronomy
}  // namespace principia
