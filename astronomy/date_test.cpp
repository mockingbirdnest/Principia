
#include "astronomy/date.hpp"

#include "gtest/gtest.h"


namespace principia {
namespace astronomy {
namespace internal_date {

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

constexpr Instant mjd0 = "1858-11-17T00:00:00Z"_TT;
constexpr Instant mjd_0 = J2000 - 51544.5 * Day;

constexpr Instant j2000_from_tt = "2000-01-01T12:00:00Z"_TT;
constexpr Instant j2000_from_tai = "2000-01-01T11:59:27,816Z"_TAI;
constexpr Instant j2000_from_utc = "2000-01-01T11:58:55,816Z"_UTC;
constexpr Instant j2000_tai = "2000-01-01T12:00:00Z"_TAI;
constexpr Instant j2000_tai_from_tt = "2000-01-01T12:00:32,184Z"_TT;

constexpr Time foo =
    "2000-01-01T11:58:55,816Z"_UTC - "2000-01-01T11:58:55,816Z"_TAI;

 constexpr Instant bat = "2000-01-01T11:58:55,816Z"_UTC;
 constexpr Instant baz = "2000-01-01T11:58:55,816Z"_TAI;

constexpr Time utc_tai =
    "2016-07-03T15:39:41Z"_UTC - "2016-07-03T15:39:41Z"_TAI;

constexpr DateTime date_Time = "1993-12-11T12:34:56,789Z"_DateTime;
constexpr DateTime date_Time_2 = "1993W496T123456.789Z"_DateTime;

constexpr TimeOfDay t = "193512,11Z"_Time;
constexpr TimeOfDay t_extended = "19:35:12,11Z"_Time;
constexpr TimeOfDay t_round = "19:35:12Z"_Time;

}  // namespace internal_date
}  // namespace astronomy
}  // namespace principia
