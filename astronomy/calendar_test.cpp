
#include "astronomy/calendar.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace astronomy {
namespace calendar {
namespace internal_calendar {

class CalendarTest : public testing::Test {};

using CalendarDeathTest = CalendarTest;

TEST_F(CalendarDeathTest, InvalidCalendarDate) {
  EXPECT_DEATH("2001-04-00"_Date, "day >= 1");
  EXPECT_DEATH("2001-02-29"_Date, "day <= month_length");
  EXPECT_DEATH("2001-03-32"_Date, "day <= month_length");
  EXPECT_DEATH("2001-04-31"_Date, "day <= month_length");
  EXPECT_DEATH("2001-00-01"_Date, "month >= 1");
  EXPECT_DEATH("2001-13-01"_Date, "month <= 12");
  EXPECT_DEATH("2001-00-01"_Date, "month >= 1");
  EXPECT_DEATH("1582-01-01"_Date, "year >= 1583");
}

TEST_F(CalendarDeathTest, InvalidTime) {
  EXPECT_DEATH("25:00:00"_Time, "hour_ <= 23");
  EXPECT_DEATH("24:01:00"_Time, "minute_ == 0");
  EXPECT_DEATH("24:00:01"_Time, "second_ == 0");
  EXPECT_DEATH("00:60:00"_Time, "minute_ <= 59");
  EXPECT_DEATH("00:00:60"_Time, "second_ <= 59");
  EXPECT_DEATH("23:59:61"_Time, "second_ == 60");
}

TEST_F(CalendarDeathTest, InvalidDateTime) {
  EXPECT_DEATH("2001-01-01T23:59:60"_DateTime,
               "date_.day.. == month_length.date_.year.., date_.month..");
}

}  // namespace internal_calendar
}  // namespace calendar
}  // namespace astronomy
}  // namespace principia
