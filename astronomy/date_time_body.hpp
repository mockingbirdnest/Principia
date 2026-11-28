#pragma once

#include "astronomy/date_time.hpp"

#include <array>
#include <iomanip>
#include <limits>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "base/macros.hpp"  // 🧙 For CONSTEXPR_CHECK.
#include "base/mod.hpp"

namespace principia {
namespace astronomy {
namespace _date_time {
namespace internal {

using namespace principia::base::_mod;

// Arithmetico-calendrical utility functions.

constexpr int mjd0_jd0_offset = 2'400'000;  // 2'400'000.5, actually.
constexpr int j2000_jd0_offset = 2'451'545;

constexpr std::array<int, 12> non_leap_year_month_lengths{
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

constexpr bool is_leap_year(int const year, Calendar const calendar) {
  return calendar == Calendar::Julian
             ? year % 4 == 0
             : year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

constexpr int year_length(int const year, Calendar const calendar) {
  return is_leap_year(year, calendar) ? 366 : 365;
}

constexpr int month_length(int const year, int const month,
                           Calendar const calendar) {
  return (is_leap_year(year, calendar) && month == 2)
             ? 29
             : non_leap_year_month_lengths[month - 1];
}

// Result in [1, 7], 1 is Monday.
constexpr int day_of_week_on_january_1st(int const year) {
  // See [Mee98, p. 65].
  return mod(Date::Calendar(year, 1, 1).jd() + 1.5, 7);
}

constexpr int number_of_iso_weeks_in_year(int const year) {
  return day_of_week_on_january_1st(year) == 4 ||
                 (is_leap_year(year, Calendar::Gregorian) &&
                  day_of_week_on_january_1st(year) == 3)
             ? 53
             : 52;
}

// Returns the ordinal in `year` of the first day of the first week of `year`.
// The result is in [-2, 4], with values in [-2, 0] meaning that the first week
// of `year` starts in `year - 1`.
// A result in [-2, 1] means that the first day of `year` is in the first week
// of `year`; otherwise, it is in the last week of `year - 1`.
constexpr int ordinal_of_w_01_1(int const year) {
  return mod(2 - day_of_week_on_january_1st(year), 7, -2);
}

// Returns the number formed by taking `end - begin` increasingly significant
// digits, starting from the digit of the (10 ** `begin`)s.
constexpr std::int64_t digit_range(std::int64_t const digits,
                                   int const begin,
                                   int const end) {
  CONSTEXPR_CHECK(digits >= 0);
  CONSTEXPR_CHECK(begin >= 0);
  CONSTEXPR_CHECK(begin <= end);
  if (begin == end) {
    return 0;
  } else if (begin == 0) {
    return digit_range(digits / 10, 0, end - 1) * 10 + digits % 10;
  } else {
    return digit_range(digits / 10, begin - 1, end - 1);
  }
}

// Returns x * 10 ** count.
constexpr std::int64_t shift_left(std::int64_t const x, int const count) {
  CONSTEXPR_CHECK(count >= 0);
  return count == 0 ? x : shift_left(x * 10, count - 1);
}

// Returns x / 10 ** count.
constexpr std::int64_t shift_right(std::int64_t const x, int const count) {
  CONSTEXPR_CHECK(count >= 0);
  return count == 0 ? x : shift_right(x / 10, count - 1);
}

// Implementation of class `Date`.

constexpr Date Date::YYYYMMDD(
    std::int64_t const digits,
    std::optional<_date_time::Calendar> const calendar) {
  auto const abs_digits = digits < 0 ? -digits : digits;
  auto const sign = digits < 0 ? -1 : 1;
  CONSTEXPR_CHECK(abs_digits <= 9999'99'99);
  return Date::Calendar(sign * digit_range(abs_digits, 4, 8),
                        digit_range(abs_digits, 2, 4),
                        digit_range(abs_digits, 0, 2),
                        calendar);
}

constexpr Date Date::YYYYDDD(
    std::int64_t const digits,
    std::optional<_date_time::Calendar> const calendar) {
  auto const sign = digits < 0 ? -1 : 1;
  CONSTEXPR_CHECK(digits <= 9999'999);
  return Date::Ordinal(sign * digit_range(digits, 3, 7),
                       digit_range(digits, 0, 3),
                       calendar);
}

constexpr Date Date::YYYYwwD(std::int64_t const digits) {
  CONSTEXPR_CHECK(digits >= 0);
  CONSTEXPR_CHECK(digits <= 9999'99'9);
  return Date::Week(digit_range(digits, 3, 7),
                    digit_range(digits, 1, 3),
                    digit_range(digits, 0, 1));
}

constexpr Date Date::Calendar(int const year, int const month, int const day,
                              std::optional<_date_time::Calendar> calendar) {
  if (!calendar.has_value()) {
    CONSTEXPR_CHECK(year >= 1583);
    calendar = _date_time::Calendar::Gregorian;
  }
  CONSTEXPR_CHECK(year <= 9999);
  CONSTEXPR_CHECK(month >= 1);
  CONSTEXPR_CHECK(month <= 12);
  CONSTEXPR_CHECK(day >= 1);
  CONSTEXPR_CHECK(day <= month_length(year, month, *calendar));
  return Date(year, month, day, *calendar);
}

constexpr Date Date::Ordinal(int const year, int const day,
                             std::optional<_date_time::Calendar> calendar) {
  if (!calendar.has_value()) {
    CONSTEXPR_CHECK(year >= 1583);
    calendar = _date_time::Calendar::Gregorian;
  }
  CONSTEXPR_CHECK(day >= 1);
  CONSTEXPR_CHECK(day <= year_length(year, *calendar));
  return Date::JD(Date(year, 1, 1, *calendar).jd() + (day - 1));
}

constexpr Date Date::Week(int const year, int const week, int const day) {
  CONSTEXPR_CHECK(year >= 1583);
  CONSTEXPR_CHECK(week >= 1);
  CONSTEXPR_CHECK(week <= number_of_iso_weeks_in_year(year));
  CONSTEXPR_CHECK(day >= 1);
  CONSTEXPR_CHECK(day <= 7);
  return Date::JD(Date::Calendar(year, 1, 1).jd() +
                         (ordinal_of_w_01_1(year) - 1) +
                         (week - 1) * 7 + (day - 1));
}

constexpr Date Date::JD(double jd) {
  // The calculation and the notation follow [Mee98, p. 63].
  // We use casting to std::int64_t as a constexpr std::trunc.  This corresponds
  // to Meeus’s INT.
  double const Z = static_cast<std::int64_t>(jd + 0.5);
  CONSTEXPR_CHECK(Z == jd + 0.5);
  // We require that jd represent midnight, so that the fractional part of jd +
  // 0.5 must be 0.
  constexpr double F = 0;
  double A = std::numeric_limits<double>::quiet_NaN();
  _date_time::Calendar calendar{};
  if (Z < 2299'161) {
    A = Z;
    calendar = Calendar::Julian;
  } else {
    double const α = static_cast<std::int64_t>((Z - 1867'216.25) / 36524.25);
    A = Z + 1 + α - static_cast<std::int64_t>(α / 4);
    calendar = Calendar::Gregorian;
  }
  double const B = A + 1524;
  double const C = static_cast<std::int64_t>((B - 122.1) / 365.25);
  double const D = static_cast<std::int64_t>(365.25 * C);
  double const E = static_cast<std::int64_t>((B - D) / 30.6001);
  int const d = B - D - static_cast<std::int64_t>(30.6001 * E) + F;
  int const m = E < 14 ? E - 1 : E - 13;
  int const y = m > 2 ? C - 4716 : C - 4715;
  return Date(y, m, d, calendar);
}

constexpr int Date::year() const {
  return year_;
}

constexpr int Date::month() const {
  return month_;
}

constexpr int Date::day() const {
  return day_;
}

constexpr _date_time::Calendar Date::calendar() const {
  return calendar_;
}

constexpr int Date::ordinal() const {
  return jd() - Date(year_, 1, 1, calendar_).jd() + 1;
}

constexpr double Date::jd() const {
  // The calculation and the notation follow [Mee98, p. 60 sq.].
  double Y = year_;
  double M = month_;
  double const D = day_;
  if (M <= 2) {
    // January and February are considered to be the 13th and 14th months of the
    // preceding year.
    Y = Y - 1;
    M = M + 12;
  }
  double B = std::numeric_limits<double>::quiet_NaN();
  if (calendar_ == Calendar::Julian) {
     B = 0;
  } else {
    double const A = static_cast<std::int64_t>(Y / 100);
    B = 2 - A + static_cast<std::int64_t>(A / 4);
  }
  return static_cast<std::int64_t>(365.25 * (Y + 4716)) +
         static_cast<std::int64_t>(30.6001 * (M + 1)) + D + B - 1524.5;
}

constexpr int Date::mjd() const {
  // See [Mee98, p. 63].
  return jd() - 2400'000.5;
}

constexpr Date Date::next_day() const {
  return Date::JD(jd() + 1);
}

constexpr Date::Date(int const year, int const month, int const day,
                     _date_time::Calendar const calendar)
      : year_(year),
        month_(month),
        day_(day),
        calendar_(calendar) {}

// Implementation of class TimeOfDay.

constexpr TimeOfDay TimeOfDay::hhmmss_ms(int const hhmmss, int ms) {
  CONSTEXPR_CHECK(hhmmss >= 0);
  CONSTEXPR_CHECK(hhmmss <= 99'99'99);
  return TimeOfDay(digit_range(hhmmss, 4, 6),
              digit_range(hhmmss, 2, 4),
              digit_range(hhmmss, 0, 2),
              ms);
}

constexpr int TimeOfDay::hour() const {
  return hour_;
}

constexpr int TimeOfDay::minute() const {
  return minute_;
}

constexpr int TimeOfDay::second() const {
  return second_;
}

constexpr int TimeOfDay::millisecond() const {
  return millisecond_;
}

constexpr bool TimeOfDay::is_leap_second() const {
  return second_ == 60;
}

constexpr bool TimeOfDay::is_end_of_day() const {
  return hour_ == 24;
}

constexpr TimeOfDay::TimeOfDay(int const hour,
                               int const minute,
                               int const second,
                               int const millisecond)
    : hour_(hour),
      minute_(minute),
      second_(second),
      millisecond_(millisecond) {
  CONSTEXPR_CHECK(
      (hour_ == 24 && minute_ == 0 && second_ == 0 && millisecond_ == 0) ||
      ((millisecond_ >= 0 && millisecond_ <= 999) &&
       ((hour_ == 23 && minute_ == 59 && second_ == 60) ||
        (hour_ >= 0 && hour_ <= 23 &&
         minute_ >= 0 && minute_ <= 59 &&
         second_ >= 0 && second_ <= 59))));
}

// Implementation of class DateTime.

constexpr DateTime DateTime::BeginningOfDay(Date const & date) {
  return DateTime(date, TimeOfDay::hhmmss_ms(00'00'00, 0));
}

constexpr DateTime DateTime::Parse(std::string_view s) {
  // Given correctness of the date and time parts of the string, this check
  // ensures that either both are in basic format or both are in extended
  // format.
  CONSTEXPR_CHECK(s.contains('-') == s.contains(':'));
  int const index_of_T = s.find_first_of('T');
  CONSTEXPR_CHECK(index_of_T != std::string_view::npos);
  return DateTime(Date::Parse(s.substr(0, index_of_T)),
                  TimeOfDay::Parse(s.substr(index_of_T + 1)));
}

constexpr Date const& DateTime::date() const {
  return date_;
}

constexpr TimeOfDay const& DateTime::time() const {
  return time_;
}

constexpr DateTime DateTime::normalized_end_of_day() const {
  return time_.is_end_of_day() ? BeginningOfDay(date_.next_day()) : *this;
}

constexpr DateTime::DateTime(Date const date, TimeOfDay const time)
    : date_(date),
      time_(time) {
  CONSTEXPR_CHECK(!time_.is_leap_second() ||
                  date_.day() == month_length(date_.year(),
                                              date_.month(),
                                              date_.calendar()));
}

// Implementation of class JulianDate.

constexpr JulianDate JulianDate::JD(std::int64_t const digits,
                                    std::int64_t const digit_count,
                                    std::int64_t const fractional_digit_count) {
  auto const day =
      digit_range(digits, fractional_digit_count, digit_count);
  auto const fraction_numerator =
      digit_range(digits, 0, fractional_digit_count);
  auto const fraction_denominator = shift_left(1, fractional_digit_count);
  return JulianDate(day - j2000_jd0_offset,
                    fraction_numerator,
                    fraction_denominator);
}

constexpr JulianDate JulianDate::MJD(
    std::int64_t const digits,
    std::int64_t const digit_count,
    std::int64_t const fractional_digit_count) {
  auto const day =
      digit_range(digits, fractional_digit_count, digit_count);
  auto const fraction_numerator =
      digit_range(digits, 0, fractional_digit_count);
  auto const fraction_denominator = shift_left(1, fractional_digit_count);
  if (fraction_denominator == 1) {
    // The only power of 10 that's not a multiple of 2.
    CONSTEXPR_CHECK(fraction_numerator == 0);
    return JulianDate(day - j2000_jd0_offset + mjd0_jd0_offset,
                      5,
                      10);
  } else if (fraction_numerator >= fraction_denominator / 2) {
    return JulianDate(day - j2000_jd0_offset + mjd0_jd0_offset + 1,
                      fraction_numerator - fraction_denominator / 2,
                      fraction_denominator);
  } else {
    return JulianDate(day - j2000_jd0_offset + mjd0_jd0_offset,
                      fraction_numerator + fraction_denominator / 2,
                      fraction_denominator);
  }
}

constexpr std::int64_t JulianDate::day() const {
  return day_;
}

constexpr std::int64_t JulianDate::fraction_numerator() const {
  return fraction_numerator_;
}

constexpr std::int64_t JulianDate::fraction_denominator() const {
  return fraction_denominator_;
}

constexpr JulianDate::JulianDate(std::int64_t day,
                                 std::int64_t fraction_numerator,
                                 std::int64_t fraction_denominator)
    : day_(day),
      fraction_numerator_(fraction_numerator),
      fraction_denominator_(fraction_denominator) {
  CONSTEXPR_CHECK(fraction_numerator >= 0);
  CONSTEXPR_CHECK(fraction_numerator < fraction_denominator);
}

// Date parsing.

// A `DateParser` contains information about a string necessary to interpret it
// as a date representation.
class DateParser final {
 public:
  // Returns a `Date` corresponding to the representation `str`.
  // Fails unless `str` is a date representation of one of the following forms:
  // [YYYY-MM-DD], [YYYYMMDD], [YYYY-Www-D], [YYYYWwwD], [YYYY-DDD], [YYYYDDD],
  // with an optional prefix [J] or [G] and an optional sign.
  static constexpr Date Parse(std::string_view s);

 private:
  constexpr DateParser(std::int64_t digits,
                       int digit_count,
                       int hyphens,
                       int first_hyphen_index,
                       int second_hyphen_index,
                       bool has_w,
                       int w_index);

  // Returns a `DateParser` describing the given string. Fails if the string
  // does not exclusively consist of:
  //   - any number of decimal digits;
  //   - at most two hyphens;
  //   - at most one 'W'.
  static constexpr DateParser ReadToEnd(std::string_view s);

  // Returns a `Date` corresponding to the string that `*this` describes.
  // Fails if the format is invalid or the string represents an invalid date.
  constexpr Date ToDate(std::optional<Calendar> calendar, bool negative) const;

  // The number formed by all digits in the string.
  std::int64_t const digits_;
  // The number of digits.
  int const digit_count_;
  // The number of hyphens.
  int const hyphens_;
  // The index of the first hyphen.
  int const first_hyphen_index_;
  // The index of the second hyphen.
  int const second_hyphen_index_;
  // Whether the string contains a W.
  bool const has_w_;
  // The index of the W.
  int const w_index_;
};

constexpr Date DateParser::Parse(std::string_view const s) {
  std::optional<Calendar> calendar;
  std::string_view tail = s;

  if (tail[0] == 'J' || tail[0] == 'G') {
    calendar = static_cast<Calendar>(tail[0]);
    tail.remove_prefix(1);
  }
  bool negative = false;
  if (tail[0] == '+' || tail[0] == '-') {
    negative = tail[0] == '-';
    tail.remove_prefix(1);
  }
  return ReadToEnd(tail).ToDate(calendar, negative);
}

constexpr DateParser::DateParser(std::int64_t const digits,
                                 int const digit_count,
                                 int const hyphens,
                                 int const first_hyphen_index,
                                 int const second_hyphen_index,
                                 bool const has_w,
                                 int const w_index)
    : digits_(digits),
      digit_count_(digit_count),
      hyphens_(hyphens),
      first_hyphen_index_(first_hyphen_index),
      second_hyphen_index_(second_hyphen_index),
      has_w_(has_w),
      w_index_(w_index) {}

constexpr DateParser DateParser::ReadToEnd(std::string_view const s) {
  std::int64_t digits = 0;
  int digit_count = 0;
  int hyphens = 0;
  int first_hyphen_index = -1;
  int second_hyphen_index = -1;
  bool has_w = false;
  int w_index = -1;

  for (auto it = s.begin(); it != s.end(); ++it) {
    switch (*it) {
      case '-':
        CONSTEXPR_CHECK(hyphens < 2);
        if (hyphens == 0) {
          first_hyphen_index = it - s.begin();
        } else {
          second_hyphen_index = it - s.begin();
        }
        ++hyphens;
        break;
      case 'W':
        CONSTEXPR_CHECK(!has_w);
        has_w = true;
        w_index = it - s.begin();
        break;
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        digits = digits * 10 + *it - '0';
        ++digit_count;
        break;
      default:
        CONSTEXPR_CHECK(false);
        return DateParser{0, 0, 0, 0, 0, false, 0};
    }
  }
  return DateParser{digits,
                    digit_count,
                    hyphens,
                    first_hyphen_index,
                    second_hyphen_index,
                    has_w,
                    w_index};
}

constexpr Date DateParser::ToDate(std::optional<Calendar> const calendar,
                                  bool const negative) const {
  auto const signed_digits = negative ? -digits_ : digits_;
  if (digit_count_ == 8) {
    CONSTEXPR_CHECK(hyphens_ == 0 ||
                    (hyphens_ == 2 &&
                     first_hyphen_index_ == 4 &&
                     second_hyphen_index_ == 7));
    return Date::YYYYMMDD(signed_digits, calendar);
  } else {
    CONSTEXPR_CHECK(digit_count_ == 7);
    if (has_w_) {
      CONSTEXPR_CHECK((hyphens_ == 0 && w_index_ == 4) ||
                      (hyphens_ == 2 &&
                       first_hyphen_index_ == 4 &&
                       w_index_ == 5 &&
                       second_hyphen_index_ == 8));
      CONSTEXPR_CHECK(!calendar.has_value());
      return Date::YYYYwwD(signed_digits);
    } else {
      CONSTEXPR_CHECK(hyphens_ == 0 ||
            (hyphens_ == 1 &&
             first_hyphen_index_ == 4));
      return Date::YYYYDDD(signed_digits, calendar);
    }
  }
}

constexpr Date Date::Parse(std::string_view const s) {
  return DateParser::Parse(s);
}

// TimeOfDay parsing.

// A `TimeParser` contains information about a string necessary to interpret it
// as a time representation.
class TimeParser final {
 public:
  // Returns a `TimeOfDay` corresponding to the representation `str`.
  // Fails unless `str` is a valid time representation of one of the following
  // forms: [hh:mm:ss], [hhmmss], [hh:mm:ss.ss̲], [hh:mm:ss,ss̲], [hhmmss.ss̲],
  // [hhmmss,ss̲], with at most three digits after the decimal mark.
  static constexpr TimeOfDay Parse(std::string_view s);

 private:
  constexpr TimeParser(std::int64_t digits,
                       int digit_count,
                       int colons,
                       int first_colon_index,
                       int second_colon_index,
                       bool has_decimal_mark,
                       int decimal_mark_index);

  // Returns a `TimeParser` describing the given string. Fails if the string
  // does not exclusively consist of:
  // Fails if the string does not exclusively consist of:
  //   - any number of decimal digits;
  //   - at most two colons;
  //   - at most one decimal mark ('.' or ',').
  static constexpr TimeParser ReadToEnd(std::string_view s);

  // Returns a `TimeOfDay` corresponding to the string that `*this` describes.
  // Fails if the format is invalid or the string represents an invalid time.
  constexpr TimeOfDay ToTime() const;

  // The number formed by all digits in the string.
  std::int64_t const digits_;
  // The number of digits.
  int const digit_count_;
  // The number of colons.
  int const colons_;
  // The index of the first colon.
  int const first_colon_index_;
  // The index of the second colon.
  int const second_colon_index_;
  // Whether the string contains a decimal mark.
  bool const has_decimal_mark_;
  // The index of the decimal mark.
  int const decimal_mark_index_;
};

constexpr TimeOfDay TimeParser::Parse(std::string_view const s) {
  return ReadToEnd(s).ToTime();
}

constexpr TimeParser::TimeParser(std::int64_t const digits,
                                 int const digit_count,
                                 int const colons,
                                 int const first_colon_index,
                                 int const second_colon_index,
                                 bool const has_decimal_mark,
                                 int const decimal_mark_index)
    : digits_(digits),
      digit_count_(digit_count),
      colons_(colons),
      first_colon_index_(first_colon_index),
      second_colon_index_(second_colon_index),
      has_decimal_mark_(has_decimal_mark),
      decimal_mark_index_(decimal_mark_index) {}

constexpr TimeParser TimeParser::ReadToEnd(std::string_view const s) {
  std::int64_t digits = 0;
  int digit_count = 0;
  int colons = 0;
  int first_colon_index = -1;
  int second_colon_index = -1;
  bool has_decimal_mark = false;
  int decimal_mark_index = -1;

  for (auto it = s.begin(); it != s.end(); ++it) {
    switch (*it) {
      case ':':
        CONSTEXPR_CHECK(colons < 2);
        if (colons == 0) {
          first_colon_index = it - s.begin();
        } else {
          second_colon_index = it - s.begin();
        }
        ++colons;
        break;
      case ',':
      case '.':
        CONSTEXPR_CHECK(!has_decimal_mark);
        has_decimal_mark = true;
        decimal_mark_index = it - s.begin();
        break;
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        digits = digits * 10 + *it - '0';
        ++digit_count;
        break;
      default:
        CONSTEXPR_CHECK(false);
        return TimeParser{0, 0, 0, 0, 0, false, 0};
    }
  }
  return TimeParser{digits,
                    digit_count,
                    colons,
                    first_colon_index,
                    second_colon_index,
                    has_decimal_mark,
                    decimal_mark_index};
}

constexpr TimeOfDay TimeParser::ToTime() const {
  // Length of the hhmmss part before the decimal mark (after stripping colons).
  constexpr int hhmmss = 6;
  CONSTEXPR_CHECK(digit_count_ >= hhmmss);
  CONSTEXPR_CHECK(digit_count_ <= 15);
  CONSTEXPR_CHECK(colons_ == 0 ||
                  (colons_ == 2 &&
                   first_colon_index_ == 2 &&
                   second_colon_index_ == 5));
  CONSTEXPR_CHECK(
      (digit_count_ == hhmmss && !has_decimal_mark_) ||
      (has_decimal_mark_ && ((colons_ == 0 && decimal_mark_index_ == 6) ||
                             (colons_ != 0 && decimal_mark_index_ == 8))));
  return TimeOfDay::hhmmss_ms(
      digit_range(digits_, digit_count_ - hhmmss, digit_count_),
      shift_left(digit_range(digits_, 0, digit_count_ - hhmmss),
                 3 - (digit_count_ - hhmmss)));
}

constexpr TimeOfDay TimeOfDay::Parse(std::string_view const s) {
  return TimeParser::Parse(s);
}

// Julian date parsing.

// A `JulianDateParser` contains information about a string necessary to
// interpret it as a Julian date.
class JulianDateParser final {
 public:
  // Returns a `JulianDate` corresponding to the representation `str`.
  // Fails unless `str` is a valid time representation of the form [ddd] or
  // [ddd.fff].

  // Returns a `JulianDate` object corresponding to the given string interpreted
  // a Julian Date or a Modified Julian Date, respectively.
  static constexpr JulianDate ParseJD(std::string_view s);
  static constexpr JulianDate ParseMJD(std::string_view s);

 private:
  constexpr JulianDateParser(std::int64_t digits,
                             int digit_count,
                             int decimal_mark_index);

  // Returns a `JulianDateParser` describing the given string. Fails if the
  // string is not of the form [ddd] or [ddd.fff] or if it has too many digits
  // to fit in a std::int64_t.
  static constexpr JulianDateParser ReadToEnd(std::string_view s);

  // Returns a `JulianDate` corresponding to the string that `*this` describes.
  // Fails if the format is invalid.
  constexpr JulianDate ToJD() const;
  constexpr JulianDate ToMJD() const;

  // The number formed by all digits in the string.
  std::int64_t const digits_;
  // The number of digits.
  int const digit_count_;
  // The index of the decimal mark.
  int const decimal_mark_index_;
};

constexpr JulianDate JulianDateParser::ParseJD(std::string_view const s) {
  return ReadToEnd(s).ToJD();
}

constexpr JulianDate JulianDateParser::ParseMJD(std::string_view const s) {
  return ReadToEnd(s).ToMJD();
}

constexpr JulianDateParser::JulianDateParser(std::int64_t const digits,
                                             int const digit_count,
                                             int const decimal_mark_index)
    : digits_(digits),
      digit_count_(digit_count),
      decimal_mark_index_(decimal_mark_index) {}

constexpr JulianDateParser JulianDateParser::ReadToEnd(
    std::string_view const s) {
  std::int64_t digits = 0;
  int digit_count = 0;
  bool has_decimal_mark = false;
  int decimal_mark_index = -1;

  for (auto it = s.begin(); it != s.end(); ++it) {
    switch (*it) {
      case '.':
        CONSTEXPR_CHECK(!has_decimal_mark);
        has_decimal_mark = true;
        decimal_mark_index = it - s.begin();
        break;
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        CONSTEXPR_CHECK(digits <=
                        std::numeric_limits<std::int64_t>::max() / 10 - 9);
        digits = digits * 10 + *it - '0';
        ++digit_count;
        break;
      default:
        CONSTEXPR_CHECK(false);
        return JulianDateParser{0, 0, 0};
    }
  }
  return JulianDateParser{digits, digit_count, decimal_mark_index};
}

constexpr JulianDate JulianDateParser::ToJD() const {
  return JulianDate::JD(
      digits_,
      digit_count_,
      decimal_mark_index_ < 0 ? 0 : digit_count_ - decimal_mark_index_);
}

constexpr JulianDate JulianDateParser::ToMJD() const {
  return JulianDate::MJD(
      digits_,
      digit_count_,
      decimal_mark_index_ < 0 ? 0 : digit_count_ - decimal_mark_index_);
}

constexpr JulianDate JulianDate::Parse(std::string_view const s) {
  static constexpr std::string_view jd = "JD";
  static constexpr std::string_view mjd = "MJD";
  if (s.starts_with(jd)) {
    return JulianDateParser::ParseJD(s.substr(jd.size()));
  } else {
    CONSTEXPR_CHECK(s.starts_with(mjd));
    return JulianDateParser::ParseMJD(s.substr(mjd.size()));
  }
}

// Operators.

constexpr bool operator==(Date const& left, Date const& right) {
  return left.jd() == right.jd();
}

constexpr bool operator!=(Date const& left, Date const& right) {
  return !(left == right);
}

constexpr auto operator<=>(Date const& left, Date const& right) {
  return left.jd() <=> right.jd();
}

consteval Date operator""_Date(char const* const str, std::size_t const size) {
  return Date::Parse({str, size});
}

inline std::ostream& operator<<(std::ostream& out, Date const& date) {
  char const fill = out.fill();
  return out << (date.calendar() == Calendar::Julian ? "J"
                 : date.year() <= 1582               ? "G"
                                                     : "")
             << (date.year() < 0 ? "-" : "") << std::setfill('0')
             << std::setw(4) << std::abs(date.year()) << "-" << std::setw(2)
             << date.month() << "-" << std::setw(2) << date.day()
             << std::setfill(fill);
}

constexpr bool operator==(TimeOfDay const& left, TimeOfDay const& right) {
  return left.hour() == right.hour() &&
         left.minute() == right.minute() &&
         left.second() == right.second() &&
         left.millisecond() == right.millisecond();
}

constexpr bool operator!=(TimeOfDay const& left, TimeOfDay const& right) {
  return !(left == right);
}

consteval TimeOfDay operator""_Time(char const* const str,
                                    std::size_t const size) {
  return TimeOfDay::Parse({str, size});
}

inline std::ostream& operator<<(std::ostream& out, TimeOfDay const& time) {
  char const fill = out.fill();
  out << std::setfill('0') << std::setw(2) << time.hour() << ":" << std::setw(2)
      << time.minute() << ":" << std::setw(2) << time.second();
  if (time.millisecond() > 0) {
    out << "," << std::setw(3) << time.millisecond();
  }
  return out << std::setfill(fill);
}

constexpr bool operator==(DateTime const& left, DateTime const& right) {
  return left.normalized_end_of_day().date() ==
             right.normalized_end_of_day().date() &&
         left.normalized_end_of_day().time() ==
             right.normalized_end_of_day().time();
}

constexpr bool operator!=(DateTime const& left, DateTime const& right) {
  return !(left == right);
}

consteval DateTime operator""_DateTime(char const* const str,
                                       std::size_t const size) {
  return DateTime::Parse({str, size});
}

inline std::ostream& operator<<(std::ostream& out, DateTime const& date_time) {
  return out << date_time.date() << "T" << date_time.time();
}

constexpr bool IsJulian(std::string_view const s) {
  return s.starts_with("JD") || s.starts_with("MJD");
}

consteval JulianDate operator""_Julian(char const* const str,
                                       std::size_t const size) {
  return JulianDate::Parse({str, size});
}

}  // namespace internal
}  // namespace _date_time
}  // namespace astronomy
}  // namespace principia
