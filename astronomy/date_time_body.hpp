
#pragma once

#include "astronomy/date_time.hpp"

#include <array>
#include <limits>

#include "base/macros.hpp"
#include "base/mod.hpp"
#include "glog/logging.h"

namespace principia {
namespace astronomy {
namespace date_time {
namespace internal_date_time {

using base::mod;

// Arithmetico-calendrical utility functions.

constexpr int mjd0_yyyy = 1858;
constexpr int mjd0_yyyymmdd = 1858'11'17;
constexpr int mjd0_jd0_offset = 2'400'000;  // 2'400'000.5, actually.
constexpr int j2000_jd0_offset = 2'451'545;

constexpr std::array<int, 12> non_leap_year_month_lengths{
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

constexpr bool is_gregorian_leap_year(int const year) {
  return (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
}

constexpr int gregorian_year_length(int const year) {
  return is_gregorian_leap_year(year) ? 366 : 365;
}

constexpr int month_length(int year, int month) {
  return (is_gregorian_leap_year(year) && month == 2)
             ? 29
             : non_leap_year_month_lengths[month - 1];
}

// Result in [1, 7], 1 is Monday.
constexpr int day_of_week_on_january_1st(int const year) {
  // Gauss's formula, see
  // https://en.wikipedia.org/wiki/Determination_of_the_day_of_the_week#Gauss.27s_algorithm.
  return mod(1 + 5 * ((year - 1) % 4) +
                 4 * ((year - 1) % 100) +
                 6 * ((year - 1) % 400),
             7,
             1);
}

constexpr int number_of_weeks_in_year(int const year) {
  return day_of_week_on_january_1st(year) == 4 ||
                 (is_gregorian_leap_year(year) &&
                  day_of_week_on_january_1st(year) == 3)
             ? 53
             : 52;
}

// Returns the ordinal in |year| of the first day of the first week of |year|.
// The result is in [-2, 4], with values in [-2, 0] meaning that the first week
// of |year| starts in |year - 1|.
// A result in [-2, 1] means that the first day of |year| is in the first week
// of |year|; otherwise, it is in the last week of |year - 1|.
constexpr int ordinal_of_w_01_1(int const year) {
  return mod(2 - day_of_week_on_january_1st(year), 7, -2);
}

constexpr int days_in_1_year = 365;
constexpr int days_in_4_years = days_in_1_year * 4 + 1;
constexpr int days_in_100_years = days_in_4_years * 25 - 1;
constexpr int days_in_400_years = days_in_100_years * 4 + 1;

// Given the number of days |d| since 0000-01-01 (proleptic Gregorian), returns
// the Gregorian year.
constexpr int gregorian_days_from_0000_01_01_to_year(int const d) {
  CONSTEXPR_CHECK(d > 0);

  // NOTE(egg): in order to extend this to the whole proleptic Gregorian
  // calendar (including d ≤ 0), we would need to use |mod| and a division
  // consistent with it.  However, in astronomy the proleptic Julian calendar is
  // used before 1582, and that is not allowed by ISO 8601, so for now let us
  // ignore the problem and assume that there are no dates before 1583-01-01.
  int const modulo_400_years = (d - 1) % days_in_400_years;
  int const modulo_100_years = modulo_400_years % days_in_100_years;
  int const modulo_4_years = modulo_100_years % days_in_4_years;

  int const multiples_of_400_years =
      ((d - 1) / days_in_400_years) * 400;
  int const multiples_of_100_years =
      (modulo_400_years / days_in_100_years) * 100;
  int const multiples_of_4_years =
      (modulo_100_years / days_in_4_years) * 4;
  int const multiples_of_1_year =
      modulo_4_years / days_in_1_year;

  return multiples_of_400_years + multiples_of_100_years +
         multiples_of_4_years + multiples_of_1_year;
}

// Given the number of days |d| since 0000-01-01 (proleptic Gregorian), returns
// the ordinal in the current Gregorian year.
constexpr int gregorian_days_from_0000_01_01_to_ordinal(int const d) {
  CONSTEXPR_CHECK(d > 0);
  int const modulo_400_years = (d - 1) % days_in_400_years;
  int const modulo_100_years = modulo_400_years % days_in_100_years;
  int const modulo_4_years = modulo_100_years % days_in_4_years;
  int const modulo_1_year = modulo_4_years % days_in_1_year;
  return modulo_1_year + 1;
}

// The number of days since 0000-01-01 on the first day of |year|, in the
// proleptic Gregorian calendar.
// |gregorian_days_from_0000_01_01_to_year| is a left inverse of this function.
constexpr int gregorian_days_from_0000_01_01_at_start_of_year(int const year) {
  CONSTEXPR_CHECK(year > 0);
  return 1 + (year) * 365 +
             (year - 1) / 4 -
             (year - 1) / 100 +
             (year - 1) / 400;
}

// Returns the number formed by taking |end - begin| increasingly significant
// digits, starting from the digit of the (10 ** |begin|)s.
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

// Returns |date| advanced by the specified number of |days|. The result must be
// in the same year.
constexpr Date add_days_within_year(Date const& date, int const days) {
  CONSTEXPR_CHECK(days >= 0);
  if (days == 0) {
    return date;
  } else {
    bool const in_same_month =
        date.day() + days <= month_length(date.year(), date.month());
    if (in_same_month) {
      return Date::Calendar(date.year(), date.month(), date.day() + days);
    } else {
      CONSTEXPR_CHECK(date.month() <= 11);
      return add_days_within_year(
          Date::Calendar(date.year(), date.month() + 1, 1),
          days - month_length(date.year(), date.month()) + date.day() - 1);
    }
  }
}

// The |day|th day of some |year|.  The resulting date need not be in |year|.
constexpr Date arbitrary_ordinal(int const year, int const day) {
  return Date::Ordinal(
      gregorian_days_from_0000_01_01_to_year(
          gregorian_days_from_0000_01_01_at_start_of_year(year) + day - 1),
      gregorian_days_from_0000_01_01_to_ordinal(
          (gregorian_days_from_0000_01_01_at_start_of_year(year) + day - 1)));
}

// Implementation of class |Date|.

constexpr Date Date::YYYYMMDD(std::int64_t const digits) {
  CONSTEXPR_CHECK(digits >= 0);
  CONSTEXPR_CHECK(digits <= 9999'99'99);
  return Date::Calendar(digit_range(digits, 4, 8),
                        digit_range(digits, 2, 4),
                        digit_range(digits, 0, 2));
}

constexpr Date Date::YYYYDDD(std::int64_t const digits) {
  CONSTEXPR_CHECK(digits >= 0);
  CONSTEXPR_CHECK(digits <= 9999'999);
  return Date::Ordinal(digit_range(digits, 3, 7), digit_range(digits, 0, 3));
}

constexpr Date Date::YYYYwwD(std::int64_t const digits) {
  CONSTEXPR_CHECK(digits >= 0);
  CONSTEXPR_CHECK(digits <= 9999'99'9);
  return Date::Week(digit_range(digits, 3, 7),
                    digit_range(digits, 1, 3),
                    digit_range(digits, 0, 1));
}

constexpr Date Date::Calendar(int const year, int const month, int const day) {
  CONSTEXPR_CHECK(year >= 1583);
  CONSTEXPR_CHECK(year <= 9999);
  CONSTEXPR_CHECK(month >= 1);
  CONSTEXPR_CHECK(month <= 12);
  CONSTEXPR_CHECK(day >= 1);
  CONSTEXPR_CHECK(day <= month_length(year, month));
  return Date(year, month, day);
}

constexpr Date Date::Ordinal(int const year, int const day) {
  CONSTEXPR_CHECK(day >= 1);
  CONSTEXPR_CHECK(gregorian_year_length(year));
  return add_days_within_year(Date::Calendar(year, 1, 1), day - 1);
}

constexpr Date Date::Week(int const year, int const week, int const day) {
  CONSTEXPR_CHECK(week >= 1);
  CONSTEXPR_CHECK(week <= number_of_weeks_in_year(year));
  CONSTEXPR_CHECK(day >= 1);
  CONSTEXPR_CHECK(day <= 7);
  return arbitrary_ordinal(year,
                           (week - 1) * 7 + day - 1 + ordinal_of_w_01_1(year));
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

constexpr int Date::ordinal() const {
  if (day_ > 1) {
    return (day_ - 1) + Date(year_, month_, 1).ordinal();
  } else if (month_ > 1) {
    return month_length(year_, month_ - 1) +
           Date(year_, month_ - 1, 1).ordinal();
  } else {
    return 1;
  }
}

constexpr int Date::mjd() const {
  return gregorian_days_from_0000_01_01_at_start_of_year(year_) + ordinal() -
         (gregorian_days_from_0000_01_01_at_start_of_year(mjd0_yyyy) +
          Date::YYYYMMDD(mjd0_yyyymmdd).ordinal());
}

constexpr Date Date::next_day() const {
  if (day_ == month_length(year_, month_)) {
    if (month_ == 12) {
      return Date(year_ + 1, 1, 1);
    } else {
      return Date(year_, month_ + 1, 1);
    }
  } else {
    return Date(year_, month_, day_ + 1);
  }
}

constexpr Date::Date(int const year,
                     int const month,
                     int const day)
      : year_(year),
        month_(month),
        day_(day) {}

// Implementation of class Time.

constexpr Time Time::hhmmss_ms(int const hhmmss, int ms) {
  CONSTEXPR_CHECK(hhmmss >= 0);
  CONSTEXPR_CHECK(hhmmss <= 99'99'99);
  return Time(digit_range(hhmmss, 4, 6),
              digit_range(hhmmss, 2, 4),
              digit_range(hhmmss, 0, 2),
              ms);
}

constexpr int Time::hour() const {
  return hour_;
}

constexpr int Time::minute() const {
  return minute_;
}

constexpr int Time::second() const {
  return second_;
}

constexpr int Time::millisecond() const {
  return millisecond_;
}

constexpr bool Time::is_leap_second() const {
  return second_ == 60;
}

constexpr bool Time::is_end_of_day() const {
  return hour_ == 24;
}

constexpr Time::Time(int const hour,
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
  return DateTime(date, Time::hhmmss_ms(00'00'00, 0));
}

constexpr Date const& DateTime::date() const {
  return date_;
}

constexpr Time const& DateTime::time() const {
  return time_;
}

constexpr DateTime DateTime::normalized_end_of_day() const {
  return time_.is_end_of_day() ? BeginningOfDay(date_.next_day()) : *this;
}

constexpr DateTime::DateTime(Date const date, Time const time)
    : date_(date),
      time_(time) {
  CONSTEXPR_CHECK(!time_.is_leap_second() ||
                  date_.day() == month_length(date_.year(), date_.month()));
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
    // The only power of 10 that's not a power of 2.
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

// Parsing utilities.

constexpr bool contains(char const* const str,
                        std::size_t const size,
                        char const c) {
  return size > 0 && (str[0] == c || contains(str + 1, size - 1, c));
}

constexpr int index_of(char const* const str,
                       std::size_t const size,
                       char const c) {
  CONSTEXPR_CHECK(size > 0);
  return str[0] == c ? 0 : (index_of(str + 1, size - 1, c) + 1);
}

constexpr bool starts_with(char const* const str,
                           std::size_t const size,
                           char const* const prefix_str,
                           std::size_t const prefix_size) {
  CONSTEXPR_CHECK(size > 0);
  CONSTEXPR_CHECK(prefix_size > 0);
  return str[0] == prefix_str[0] &&
         (prefix_size == 1 ||
          (size > 1 &&
           starts_with(str + 1, size - 1, prefix_str + 1, prefix_size - 1)));
}

class CStringIterator final {
 public:
  constexpr CStringIterator(char const* str, std::size_t size);

  constexpr bool at_end() const;
  constexpr CStringIterator next() const;
  constexpr int index() const;
  constexpr char const& operator*() const;

 private:
  constexpr CStringIterator(char const* str, char const* end, char const* it);

  char const* const str_;
  char const* const end_;
  char const* const it_;
};

constexpr CStringIterator::CStringIterator(char const* const str,
                                           std::size_t const size)
    : str_(str), end_(str + size), it_(str) {}

constexpr bool CStringIterator::at_end() const {
  return it_ == end_;
}

constexpr CStringIterator CStringIterator::next() const {
  CONSTEXPR_CHECK(!at_end());
  return CStringIterator(str_, end_, it_ + 1);
}

constexpr int CStringIterator::index() const {
  return it_ - str_;
}

constexpr char const& CStringIterator::operator*() const {
  CONSTEXPR_CHECK(!at_end());
  return *it_;
}

constexpr CStringIterator::CStringIterator(char const* const str,
                                           char const* const end,
                                           char const* const it)
    : str_(str),
      end_(end),
      it_(it) {}

// Date parsing.

// A |DateParser| contains information about a string necessary to interpret it
// as a date representation.
class DateParser final {
 public:
  // Returns a |Date| corresponding to the representation |str|.
  // Fails unless |str| is a date representation of one of the following forms:
  // [YYYY-MM-DD], [YYYYMMDD], [YYYY-Www-D], [YYYYWwwD], [YYYY-DDD], [YYYYDDD].
  static constexpr Date Parse(char const* str, std::size_t size);

 private:
  constexpr DateParser(std::int64_t digits,
                       int digit_count,
                       int hyphens,
                       int first_hyphen_index,
                       int second_hyphen_index,
                       bool has_w,
                       int w_index);

  // Returns a |DateParser| describing the given string. Fails if the string
  // does not exclusively consist of:
  //   - any number of decimal digits;
  //   - at most two hyphens;
  //   - at most one 'W'.
  static constexpr DateParser ReadToEnd(char const* str,
                                        std::size_t size);
  static constexpr DateParser ReadToEnd(CStringIterator str,
                                        std::int64_t digits,
                                        int digit_count,
                                        int colons,
                                        int first_hyphen_index,
                                        int second_hyphen_index,
                                        bool has_w,
                                        int w_index);

  // Returns a |Date| corresponding to the string that |*this| describes.
  // Fails if the format is invalid or the string represents an invalid date.
  constexpr Date ToDate() const;

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

constexpr Date DateParser::Parse(char const* const str,
                                 std::size_t const size) {
  return ReadToEnd(str, size).ToDate();
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

constexpr DateParser DateParser::ReadToEnd(char const* const str,
                                           std::size_t const size) {
  return ReadToEnd(CStringIterator(str, size),
                   /*digits=*/0,
                   /*digit_count=*/0,
                   /*hyphens=*/0,
                   /*first_hyphen_index=*/-1,
                   /*second_hyphen_index=*/-1,
                   /*has_w=*/false,
                   /*w_index=*/-1);
}

constexpr DateParser DateParser::ReadToEnd(CStringIterator const str,
                                           std::int64_t const digits,
                                           int const digit_count,
                                           int const hyphens,
                                           int const first_hyphen_index,
                                           int const second_hyphen_index,
                                           bool const has_w,
                                           int const w_index) {
  if (str.at_end()) {
    return DateParser{digits,
                      digit_count,
                      hyphens,
                      first_hyphen_index,
                      second_hyphen_index,
                      has_w,
                      w_index};
  } else {
    switch (*str) {
      case '-':
        CONSTEXPR_CHECK(hyphens < 2);
        if (hyphens == 0) {
          return ReadToEnd(str.next(),
                           digits,
                           digit_count,
                           hyphens + 1,
                           /*first_hyphen_index=*/str.index(),
                           second_hyphen_index,
                           has_w,
                           w_index);
        } else {
          return ReadToEnd(str.next(),
                           digits,
                           digit_count,
                           hyphens + 1,
                           first_hyphen_index,
                           /*second_hyphen_index=*/str.index(),
                           has_w,
                           w_index);
        }
      case 'W':
        CONSTEXPR_CHECK(!has_w);
        return ReadToEnd(str.next(),
                         digits,
                         digit_count,
                         hyphens,
                         first_hyphen_index,
                         second_hyphen_index,
                         /*has_w=*/true,
                         /*w_index=*/str.index());
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        return ReadToEnd(str.next(),
                         digits * 10 + *str - '0',
                         digit_count + 1,
                         hyphens,
                         first_hyphen_index,
                         second_hyphen_index,
                         has_w,
                         w_index);
      default:
        CONSTEXPR_CHECK(false);
        return DateParser{0, 0, 0, 0, 0, false, 0};
    }
  }
}

constexpr Date DateParser::ToDate() const {
  if (digit_count_ == 8) {
    CONSTEXPR_CHECK(hyphens_ == 0 ||
                    (hyphens_ == 2 &&
                     first_hyphen_index_ == 4 &&
                     second_hyphen_index_ == 7));
    return Date::YYYYMMDD(digits_);
  } else {
    CONSTEXPR_CHECK(digit_count_ == 7);
    if (has_w_) {
      CONSTEXPR_CHECK((hyphens_ == 0 && w_index_ == 4) ||
                      (hyphens_ == 2 &&
                       first_hyphen_index_ == 4 &&
                       w_index_ == 5 &&
                       second_hyphen_index_ == 8));
      return Date::YYYYwwD(digits_);
    } else {
      CONSTEXPR_CHECK(hyphens_ == 0 ||
            (hyphens_ == 1 &&
             first_hyphen_index_ == 4));
      return Date::YYYYDDD(digits_);
    }
  }
}

// Time parsing.

// A |TimeParser| contains information about a string necessary to interpret it
// as a time representation.
class TimeParser final {
 public:
  // Returns a |Time| corresponding to the representation |str|.
  // Fails unless |str| is a valid time representation of one of the following
  // forms: [hh:mm:ss], [hhmmss], [hh:mm:ss.ss̲], [hh:mm:ss,ss̲], [hhmmss.ss̲],
  // [hhmmss,ss̲], with at most three digits after the decimal mark.
  static constexpr Time Parse(char const* str, std::size_t size);

 private:
  constexpr TimeParser(std::int64_t digits,
                       int digit_count,
                       int colons,
                       int first_colon_index,
                       int second_colon_index,
                       bool has_decimal_mark,
                       int decimal_mark_index);

  // Returns a |TimeParser| describing the given string. Fails if the string
  // does not exclusively consist of:
  // Fails if the string does not exclusively consist of:
  //   - any number of decimal digits;
  //   - at most two colons;
  //   - at most one decimal mark ('.' or ',').
  static constexpr TimeParser ReadToEnd(char const* str,
                                        std::size_t size);
  static constexpr TimeParser ReadToEnd(CStringIterator str,
                                        std::int64_t digits,
                                        int digit_count,
                                        int colons,
                                        int first_colon_index,
                                        int second_colon_index,
                                        bool has_decimal_mark,
                                        int decimal_mark_index);

  // Returns a |Time| corresponding to the string that |*this| describes.
  // Fails if the format is invalid or the string represents an invalid time.
  constexpr Time ToTime() const;

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

constexpr Time TimeParser::Parse(char const* const str,
                                 std::size_t const size) {
  return ReadToEnd(str, size).ToTime();
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

constexpr TimeParser TimeParser::ReadToEnd(char const* const str,
                                           std::size_t const size) {
  return ReadToEnd(CStringIterator(str, size),
                   /*digits=*/0,
                   /*digit_count*/ 0,
                   /*colons=*/0,
                   /*first_colon_index=*/-1,
                   /*second_colon_index=*/-1,
                   /*has_decimal_mark=*/false,
                   /*decimal_mark_index=*/-1);
}

constexpr TimeParser TimeParser::ReadToEnd(CStringIterator const str,
                                           std::int64_t const digits,
                                           int const digit_count,
                                           int const colons,
                                           int const first_colon_index,
                                           int const second_colon_index,
                                           bool const has_decimal_mark,
                                           int const decimal_mark_index) {
  if (str.at_end()) {
    return TimeParser(digits,
                      digit_count,
                      colons,
                      first_colon_index,
                      second_colon_index,
                      has_decimal_mark,
                      decimal_mark_index);
  } else {
    switch (*str) {
      case ':':
        CONSTEXPR_CHECK(colons < 2);
        if (colons == 0) {
          return ReadToEnd(str.next(),
                           digits,
                           digit_count,
                           colons + 1,
                           /*first_colon_index=*/str.index(),
                           second_colon_index,
                           has_decimal_mark,
                           decimal_mark_index);
        } else {
          return ReadToEnd(str.next(),
                           digits,
                           digit_count,
                           colons + 1,
                           first_colon_index,
                           /*second_colon_index=*/str.index(),
                           has_decimal_mark,
                           decimal_mark_index);
        }
      case ',':
      case '.':
        CONSTEXPR_CHECK(!has_decimal_mark);
        return ReadToEnd(str.next(),
                         digits,
                         digit_count,
                         colons,
                         first_colon_index,
                         second_colon_index,
                         /*has_decimal_mark=*/true,
                         /*decimal_mark_index=*/str.index());
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        return ReadToEnd(str.next(),
                         digits * 10 + *str - '0',
                         digit_count + 1,
                         colons,
                         first_colon_index,
                         second_colon_index,
                         has_decimal_mark,
                         decimal_mark_index);
      default:
        CONSTEXPR_CHECK(false);
        return TimeParser{0, 0, 0, 0, 0, false, 0};
    }
  }
}

constexpr Time TimeParser::ToTime() const {
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
  return Time::hhmmss_ms(
      digit_range(digits_, digit_count_ - hhmmss, digit_count_),
      shift_left(digit_range(digits_, 0, digit_count_ - hhmmss),
                 3 - (digit_count_ - hhmmss)));
}

// Julian date parsing.

// A |JulianDateParser| contains information about a string necessary to
// interpret it as a Julian date.
class JulianDateParser final {
 public:
  // Returns a |JulianDate| corresponding to the representation |str|.
  // Fails unless |str| is a valid time representation of the form [ddd] or
  // [ddd.fff].

  // Returns a |JulianDate| object corresponding to the given string interpreted
  // a Julian Date or a Modified Julian Date, respectively.
  static constexpr JulianDate ParseJD(char const* str, std::size_t size);
  static constexpr JulianDate ParseMJD(char const* str, std::size_t size);

 private:
  constexpr JulianDateParser(std::int64_t digits,
                             int digit_count,
                             bool has_decimal_mark,
                             int decimal_mark_index);

  // Returns a |JulianDateParser| describing the given string. Fails if the
  // string is not of the form [ddd] or [ddd.fff] or if it has too many digits
  // to fit in a std::int64_t.
  static constexpr JulianDateParser ReadToEnd(char const* str,
                                              std::size_t size);
  static constexpr JulianDateParser ReadToEnd(CStringIterator str,
                                              std::int64_t digits,
                                              int digit_count,
                                              bool has_decimal_mark,
                                              int decimal_mark_index);

  // Returns a |JulianDate| corresponding to the string that |*this| describes.
  // Fails if the format is invalid.
  constexpr JulianDate ToJD() const;
  constexpr JulianDate ToMJD() const;

  // The number formed by all digits in the string.
  std::int64_t const digits_;
  // The number of digits.
  int const digit_count_;
  // Whether the string contains a decimal mark.
  bool const has_decimal_mark_;
  // The index of the decimal mark.
  int const decimal_mark_index_;
};

constexpr JulianDate JulianDateParser::ParseJD(char const* const str,
                                               std::size_t const size) {
  return ReadToEnd(str, size).ToJD();
}

constexpr JulianDate JulianDateParser::ParseMJD(char const* const str,
                                                std::size_t const size) {
  return ReadToEnd(str, size).ToMJD();
}

constexpr JulianDateParser::JulianDateParser(std::int64_t const digits,
                                             int const digit_count,
                                             bool const has_decimal_mark,
                                             int const decimal_mark_index)
    : digits_(digits),
      digit_count_(digit_count),
      has_decimal_mark_(has_decimal_mark),
      decimal_mark_index_(decimal_mark_index) {}

constexpr JulianDateParser JulianDateParser::ReadToEnd(char const* const str,
                                                       std::size_t const size) {
  return ReadToEnd(CStringIterator(str, size),
                   /*digits=*/0,
                   /*digit_count*/ 0,
                   /*has_decimal_mark=*/false,
                   /*decimal_mark_index=*/-1);
}

constexpr JulianDateParser JulianDateParser::ReadToEnd(
    CStringIterator const str,
    std::int64_t const digits,
    int const digit_count,
    bool const has_decimal_mark,
    int const decimal_mark_index) {
  if (str.at_end()) {
    return JulianDateParser(digits,
                            digit_count,
                            has_decimal_mark,
                            decimal_mark_index);
  } else {
    switch (*str) {
      case '.':
        CONSTEXPR_CHECK(!has_decimal_mark);
        return ReadToEnd(str.next(),
                         digits,
                         digit_count,
                         /*has_decimal_mark=*/true,
                         /*decimal_mark_index=*/str.index());
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        CONSTEXPR_CHECK(digits <=
                        std::numeric_limits<std::int64_t>::max() / 10 - 9);
        return ReadToEnd(str.next(),
                         digits * 10 + *str - '0',
                         digit_count + 1,
                         has_decimal_mark,
                         decimal_mark_index);
      default:
        CONSTEXPR_CHECK(false);
        return JulianDateParser{0, 0, false, 0};
    }
  }
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

// Operators.

constexpr bool operator==(Date const& left, Date const& right) {
  return left.year() == right.year() &&
         left.month() == right.month() &&
         left.day() == right.day();
}

constexpr bool operator!=(Date const& left, Date const& right) {
  return !(left == right);
}

constexpr bool operator<(Date const& left, Date const& right) {
  return left.year() < right.year() ||
         (left.year() == right.year() &&
          (left.month() < right.month() ||
           (left.month() == right.month() && left.day() < right.day())));
}

constexpr bool operator>(Date const& left, Date const& right) {
  return right < left;
}

constexpr bool operator<=(Date const& left, Date const& right) {
  return !(right < left);
}

constexpr bool operator>=(Date const& left, Date const& right) {
  return !(left < right);
}

constexpr Date operator""_Date(char const* const str, std::size_t const size) {
  return DateParser::Parse(str, size);
}

constexpr bool operator==(Time const& left, Time const& right) {
  return left.hour() == right.hour() &&
         left.minute() == right.minute() &&
         left.second() == right.second() &&
         left.millisecond() == right.millisecond();
}

constexpr bool operator!=(Time const& left, Time const& right) {
  return !(left == right);
}

constexpr Time operator""_Time(char const* const str, std::size_t const size) {
  return TimeParser::Parse(str, size);
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

constexpr DateTime operator""_DateTime(char const* const str,
                                       std::size_t const size) {
  // Given correctness of the date and time parts of the string, this check
  // ensures that either both are in basic format or both are in extended
  // format.
  CONSTEXPR_CHECK(contains(str, size, '-') == contains(str, size, ':'));
  const int index_of_T = index_of(str, size, 'T');
  return DateTime(DateParser::Parse(str, index_of_T),
                  TimeParser::Parse(str + index_of_T + 1,
                                    size - (index_of_T + 1)));
}

constexpr bool IsJulian(char const* const str, std::size_t const size) {
  return starts_with(str, size, "JD", 2) || starts_with(str, size, "MJD", 3);
}

constexpr JulianDate operator""_Julian(char const* const str,
                                       std::size_t const size) {
  if (starts_with(str, size, "JD", 2)) {
    return JulianDateParser::ParseJD(str + 2, size - 2);
  } else {
    CONSTEXPR_CHECK(starts_with(str, size, "MJD", 3));
    return JulianDateParser::ParseMJD(str + 3, size - 3);
  }
}

}  // namespace internal_date_time
}  // namespace date_time
}  // namespace astronomy
}  // namespace principia
