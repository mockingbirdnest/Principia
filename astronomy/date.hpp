
#pragma once

#include <array>
#include <cstdint>

#include "glog/logging.h"

#define CHECKING(condition, expression) \
  ((condition) ? (expression) : (CHECK(condition), (expression)))

namespace principia {
namespace astronomy {
namespace internal_date {

class Date {
 public:
  static constexpr Date YYYYMMDD(std::int64_t digits);
  static constexpr Date YYYYwwD(std::int64_t digits);
  static constexpr Date YYYYDDD(std::int64_t digits);

  constexpr int year() const;
  constexpr int month() const;
  constexpr int day() const;

 private:
  constexpr Date(int const year,
                 std::int8_t const month,
                 std::int8_t const day);

  constexpr Date const& checked() const;

  int const year_;
  int const month_;
  int const day_;

  friend constexpr Date add_days(Date const& date, int const days);
  friend struct OrdinalDate;
  friend constexpr Date operator""_date(char const* string, std::size_t size);
};

// The following is horrendous because written in C++11, thus functionally; it
// should perhaps be reworked once MSVC offers proper C++14 constexpr.

constexpr Date::Date(int const year,
                     std::int8_t const month,
                     std::int8_t const day)
      : year_(year),
        month_(month),
        day_(day) {}

constexpr int Date::year() const {
  return year_;
}

constexpr int Date::month() const {
  return month_;
}

constexpr int Date::day() const {
  return day_;
}

constexpr std::int64_t digit_range(std::int64_t const digits,
                                   int const begin,
                                   int const size) {
  return CHECKING(
      digits >= 0 && begin >= 0 && size >= 0,
      (size == 0) ? 0 : ((begin == 0)
                             ? digit_range(digits / 10, 0, size - 1) * 10 +
                                   digits % 10
                             : digit_range(digits / 10, begin - 1, size)));
}

constexpr bool is_gregorian_leap_year(int const year) {
  return (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
}

constexpr int gregorian_year_length(int const year) {
  return is_gregorian_leap_year(year) ? 366 : 365;
}

constexpr std::array<int, 12>
    non_leap_year_month_lengths{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

constexpr int month_length(int year, int month) {
  return (is_gregorian_leap_year(year) && month == 2)
             ? 29
             : non_leap_year_month_lengths[month - 1];
}

constexpr int mod_7_offset_1 (int const x) {
  return x % 7 == 0 ? 7 : (x % 7);
}

// Result in [1, 7], 1 is Monday.
constexpr int day_of_week_on_january_1st(int const year) {
  // Gauss's formula, see
  // https://en.wikipedia.org/wiki/Determination_of_the_day_of_the_week#Gauss.27s_algorithm.
  return mod_7_offset_1(1 +
                        5 * ((year - 1) % 4) +
                        4 * ((year - 1) % 100) +
                        6 * ((year - 1) % 400));
}

constexpr Date Date::YYYYMMDD(std::int64_t digits) {
  return CHECKING(digits < 9999'99'99,
                  Date(digit_range(digits, 4, 4),
                       digit_range(digits, 2, 2),
                       digit_range(digits, 0, 2)).checked());
}

constexpr Date const& Date::checked() const {
  return CHECKING(year_ >= 1583 && year_ <= 9999 &&
                  month_ >= 1 && month_ <= 12 &&
                  day_ >= 1 && day_ <= month_length(year_, month_),
                  *this);
}

constexpr Date d = Date::YYYYMMDD(1993'12'11);
//constexpr Date d1 = Date::YYYYMMDD(1993'02'30);

// Returns |date| advanced by the specified number of |days|. |days| must be
// positive, and the result must be in the same year.
constexpr Date add_days(Date const& date, int const days) {
  return CHECKING(
      days > 0,
      days == 0
          ? date
          : (date.day() + days > month_length(date.year(), date.month())
                 ? CHECKING(
                       date.month() <= 11,
                       add_days(Date(date.year(), date.month() + 1, 1),
                                days - month_length(date.year(), date.month()) +
                                    date.day() - 1))
                 : Date(date.year(), date.month(), date.day() + days)));
}

struct OrdinalDate {
  constexpr OrdinalDate(int const year, int const day) : year(year), day(day) {}

  constexpr OrdinalDate const& checked() const {
    return CHECKING(year >= 1583 && year < 9999 &&
                    day >= 1 && day <= gregorian_year_length(year),
                    *this);
  }

  constexpr OrdinalDate const& normalized() const {
    return day < 1
        ? OrdinalDate(year - 1,
                      gregorian_year_length(year - 1) + day).normalized()
        : (day > gregorian_year_length(year)
               ? OrdinalDate(year + 1,
                             day - gregorian_year_length(year)).normalized()
               : *this);
  }

  constexpr Date ToDate() const {
    return add_days(Date(year, 1, 1), day - 1);
  }

  int year;
  int day;
};

constexpr Date Date::YYYYDDD(std::int64_t digits) {
  return CHECKING(
      digits > 0 && digits < 9999'999,
      OrdinalDate(digit_range(digits, 3, 4),
                  digit_range(digits, 0, 3)).checked().ToDate());
}

constexpr bool is_53_week_year(int const year) {
  return day_of_week_on_january_1st(year) == 4 ||
         (is_gregorian_leap_year(year) &&
          day_of_week_on_january_1st(year) == 3);
}

constexpr int ordinal_of_w_01_1(int const year) {
  return day_of_week_on_january_1st(year) <= 4
             ? 2 - day_of_week_on_january_1st(year)
             : (9 - day_of_week_on_january_1st(year));
}

struct WeekDate {
  constexpr WeekDate(int year, int week, int day)
      : year(year), week(week), day(day) {}

  constexpr WeekDate const& checked() const {
    return CHECKING(year >= 1583 && year < 9999 &&
                    week >= 1 && week <= (is_53_week_year(year) ? 53 : 52) &&
                    day >= 1 && day <= 7,
                    *this);
  }

  constexpr Date ToDate() const {
    return OrdinalDate(
        year,
        (week - 1) * 7 + day - 1 + ordinal_of_w_01_1(year)).
            normalized().checked().ToDate();
  }

  int year;
  int week;
  int day;
};

constexpr Date Date::YYYYwwD(std::int64_t digits) {
  return CHECKING(digits > 0 && digits < 9999'99'9,
                  WeekDate(digit_range(digits, 3, 4),
                           digit_range(digits, 1, 2),
                           digit_range(digits, 0, 1)).ToDate());
}

namespace from_integers_test {
constexpr Date egg_month = Date::YYYYMMDD(1993'12'11);
static_assert(egg_month.year() == 1993, "bad year");
static_assert(egg_month.month() == 12, "bad month");
static_assert(egg_month.day() == 11, "bad day");
constexpr Date egg_week = Date::YYYYwwD(1993'49'6);
static_assert(egg_week.year() == 1993, "bad year");
static_assert(egg_week.month() == 12, "bad month");
static_assert(egg_week.day() == 11, "bad day");
constexpr Date egg_ordinal = Date::YYYYDDD(1993'345);
static_assert(egg_ordinal.year() == 1993, "bad year");
static_assert(egg_ordinal.month() == 12, "bad month");
static_assert(egg_ordinal.day() == 11, "bad day");
}

struct DateStringInfo {
  constexpr DateStringInfo Fill() const {
    return read == size
               ? *this
               : string[read] == '-'
                     ? CHECKING(
                           hyphens < 2,
                           hyphens == 0
                               ? DateStringInfo({string,
                                                 size,
                                                 read + 1,
                                                 digits,
                                                 digit_count,
                                                 hyphens + 1,
                                                 /*first_hyphen_index=*/read,
                                                 second_hyphen_index,
                                                 has_w,
                                                 w_index}).Fill()
                               : DateStringInfo({string,
                                                 size,
                                                 read + 1,
                                                 digits,
                                                 digit_count,
                                                 hyphens + 1,
                                                 first_hyphen_index,
                                                 /*second_hyphen_index=*/read,
                                                 has_w,
                                                 w_index})).Fill()
                     : string[read] == 'W'
                           ? CHECKING(!has_w,
                                      DateStringInfo({string,
                                                      size,
                                                      read + 1,
                                                      digits,
                                                      digit_count,
                                                      hyphens,
                                                      first_hyphen_index,
                                                      second_hyphen_index,
                                                      /*has_w=*/true,
                                                      /*w_index=*/read})).Fill()
                           : CHECKING(
                                 string[read] >= '0' && string[read] <= '9',
                                 DateStringInfo(
                                     {string,
                                      size,
                                      read + 1,
                                      digits * 10 + string[read] - '0',
                                      digit_count + 1,
                                      hyphens,
                                      first_hyphen_index,
                                      second_hyphen_index,
                                      has_w,
                                      w_index}).Fill());
  }

  constexpr Date ToDate() const {
    return digit_count == 8
               ? CHECKING(
                     hyphens == 0 || (hyphens == 2 && first_hyphen_index == 4 &&
                                      second_hyphen_index == 7),
                     Date::YYYYMMDD(digits))
               : CHECKING(
                     digit_count == 7,
                     has_w ? CHECKING(
                                 (hyphens == 0 && w_index == 4) ||
                                     (hyphens == 2 && first_hyphen_index == 4 &&
                                      w_index == 5 && second_hyphen_index == 8),
                                 Date::YYYYwwD(digits))
                           : CHECKING(hyphens == 0 || (hyphens == 1 &&
                                                       first_hyphen_index == 4),
                                      Date::YYYYDDD(digits)));
  }

  char const* const string;
  std::size_t size;
  int const read;
  std::int64_t digits;
  int digit_count;
  int const hyphens;
  int const first_hyphen_index;
  int const second_hyphen_index;
  bool const has_w;
  int const w_index;
};

constexpr Date operator""_date(char const* string, std::size_t size) {
  return DateStringInfo{string,
                        size,
                        /*read=*/0,
                        /*digits=*/0,
                        /*digit_count=*/0,
                        /*hyphens=*/0,
                        /*first_hyphen_index=*/-1,
                        /*second_hyphen_index=*/-1,
                        /*has_w=*/false,
                        /*w_index=*/-1}.Fill().ToDate();
}

class TimeOfDay {
 public:
  static constexpr TimeOfDay hhmmss_ns(int const hhmmss, int ns) {
    return TimeOfDay(digit_range(hhmmss, 4, 2),
                     digit_range(hhmmss, 2, 2),
                     digit_range(hhmmss, 0, 2),
                     ns).checked();
  }

  constexpr int hour() const;
  constexpr int minute() const;
  constexpr int second() const;
  constexpr int nanosecond() const;

  constexpr bool is_leap_second() const;

 private:
  constexpr TimeOfDay(int const hour,
                      int const minute,
                      int const second,
                      int const nanosecond)
      : hour_(hour),
        minute_(minute),
        second_(second),
        nanosecond_(nanosecond) {}

  // Checks that this represents a valid time of day as per ISO 8601, thus
  // that the components are in the normal range, or that the object represents
  // a time in a leap second, or that it represents the end of the day.
  constexpr TimeOfDay const& checked() const {
    return CHECKING(
        (hour_ == 24 && minute_ == 0 && second_ == 0 && nanosecond_ == 0) ||
            ((nanosecond_ >= 0 && nanosecond_ <= 999'999'999) &&
             ((hour_ == 23 && minute_ == 59 && second_ == 60) ||
              (hour_ >= 0 && hour_ <= 23 && minute_ >= 0 && minute_ <= 59 &&
               second_ >= 0 && second_ <= 59))),
        *this);
  }

  int const hour_;
  int const minute_;
  int const second_;
  int const nanosecond_;
};

constexpr int TimeOfDay::hour() const {
  return hour_;
}

constexpr int TimeOfDay::minute() const {
  return minute_;
}

constexpr int TimeOfDay::second() const {
  return second_;
}

constexpr int TimeOfDay::nanosecond() const {
  return nanosecond_;
}

constexpr bool TimeOfDay::is_leap_second() const {
  return second_ == 60;
}

constexpr std::int64_t add_0s(std::int64_t const x, int const count) {
  return count == 0 ? x : add_0s(x * 10, count - 1);
}

struct TimeStringInfo {
  constexpr TimeStringInfo Fill() const {
    return read == size
               ? *this
               : string[read] == ':'
                     ? CHECKING(
                           colons < 2,
                           colons == 0
                               ? TimeStringInfo({string,
                                                 size,
                                                 read + 1,
                                                 digits,
                                                 digit_count,
                                                 colons + 1,
                                                 /*first_colon_index=*/read,
                                                 second_colon_index,
                                                 has_decimal_mark,
                                                 decimal_mark_index}).Fill()
                               : TimeStringInfo({string,
                                                 size,
                                                 read + 1,
                                                 digits,
                                                 digit_count,
                                                 colons + 1,
                                                 first_colon_index,
                                                 /*second_colon_index=*/read,
                                                 has_decimal_mark,
                                                 decimal_mark_index})).Fill()
                     : string[read] == ',' || string[read] == '.'
                           ? CHECKING(
                                 !has_decimal_mark,
                                 TimeStringInfo(
                                    {string,
                                     size,
                                     read + 1,
                                     digits,
                                     digit_count,
                                     colons,
                                     first_colon_index,
                                     second_colon_index,
                                     /*has_decimal_mark=*/true,
                                     /*decimal_mark_index=*/read}).Fill())
                           : string[read] == 'Z'
                                 ? CHECKING(
                                       read == size - 1,
                                       TimeStringInfo(
                                           {string,
                                            size,
                                            read + 1,
                                            digits,
                                            digit_count,
                                            colons,
                                            first_colon_index,
                                            second_colon_index,
                                            has_decimal_mark,
                                            decimal_mark_index}).Fill())
                                 : CHECKING(
                                       string[read] >= '0' &&
                                       string[read] <= '9',
                                       TimeStringInfo(
                                           {string,
                                            size,
                                            read + 1,
                                            digits * 10 + string[read] - '0',
                                            digit_count + 1,
                                            colons,
                                            first_colon_index,
                                            second_colon_index,
                                            has_decimal_mark,
                                            decimal_mark_index}).Fill());
  }

  constexpr TimeOfDay ToTime() const {
    return CHECKING(
        digit_count >= 6 &&
            (colons == 0 || (colons == 2 && first_colon_index == 2 &&
                             second_colon_index == 5)) &&
            ((digit_count == 6 && !has_decimal_mark) ||
             (has_decimal_mark &&
              ((colons == 0 && decimal_mark_index == 6) ||
               (colons != 0 && decimal_mark_index == 8)))) &&
            digit_count <= 15 &&
            string[size - 1] == 'Z',
        TimeOfDay::hhmmss_ns(digit_range(digits, digit_count - 6, 6),
                           add_0s(digit_range(digits, 0, digit_count - 6),
                                  9 - (digit_count - 6))));
  }

  char const* const string;
  std::size_t size;
  int const read;
  std::int64_t digits;
  int digit_count;
  int const colons;
  int const first_colon_index;
  int const second_colon_index;
  bool const has_decimal_mark;
  int const decimal_mark_index;
};

constexpr TimeOfDay operator""_time(char const* string, std::size_t size) {
  return TimeStringInfo{string,
                        size,
                        /*read=*/0,
                        /*digits=*/0,
                        /*digit_count=*/0,
                        /*colons=*/0,
                        /*first_colon_index=*/-1,
                        /*second_colon_index=*/-1,
                        /*has_decimal_mark=*/false,
                        /*decimal_mark_index=*/-1}.Fill().ToTime();
}

struct DateTime {
  // We use an explicit constructor here because too many initializer-lists in
  // constexpr code confuse the compiler and eventually make it crash and burn.
  // It removes one layer of brackets anyway (we're in macros), and this struct
  // has few members.
  constexpr DateTime(Date const date, TimeOfDay const time)
      : date(date),
        time(time) {}

  // Checks that |time| does not represent a leap second unless |date| is the
  // last day of June, December, March, or September.
  constexpr DateTime const& checked() const {
    return CHECKING(
        !time.is_leap_second() ||
            (date.day() == month_length(date.year(), date.month()) &&
             (date.month() == 6 || date.month() == 12 || date.month() == 3 ||
              date.month() == 9)),
        *this);
  }

  Date const date;
  TimeOfDay const time;
};

constexpr bool contains(char const* string, std::size_t size, char const c) {
  return (size > 0 && string[0] == c) || contains(string + 1, size - 1, c);
}

constexpr int index_of(char const* string, std::size_t size, char const c) {
  return CHECKING(size > 0,
                  string[0] == c ? 0 : (index_of(string + 1, size - 1, c) + 1));
}

constexpr DateTime operator""_date_time(char const* string, std::size_t size) {
  return CHECKING(
      contains(string, size, '-') == contains(string, size, ':'),
      DateTime(
          operator""_date(string, index_of(string, size, 'T')),
          operator""_time(string + index_of(string, size, 'T') + 1,
                          size - (index_of(string, size, 'T') + 1))).checked());
}

constexpr DateTime date_time = "1993-12-11T12:34:56,789Z"_date_time;

constexpr TimeOfDay t = "193512,11Z"_time;
constexpr TimeOfDay t_extended = "19:35:12,11Z"_time;
constexpr TimeOfDay t_round = "19:35:12Z"_time;

namespace basic_format_test {
constexpr Date egg_month = "19931211"_date;
constexpr Date egg_week = "1993W496"_date;
constexpr Date egg_ordinal = "1993345"_date;
}

namespace extended_format_test {
constexpr Date egg_month = "1993-12-11"_date;
constexpr Date egg_week = "1993-W49-6"_date;
constexpr Date egg_ordinal = "1993-345"_date;
}

}  // namespace internal_date
}  // namespace astronomy
}  // namespace principia
