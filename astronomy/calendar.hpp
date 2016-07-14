
#pragma once

#include <cstddef>
#include <cstdint>

#include "quantities/quantities.hpp"

// This file has an additional layer of namespacing for several reasons:
//   - we are not entirely sure whether we really want to expose this API;
//   - |astronomy::Time| for a class that represents the time of day is weird;
//   - until we use internal namespace properly everywhere, |astronomy::Time|
//     would clash with |quantities::Time| in some parts of |astronomy|.

namespace principia {
namespace astronomy {
namespace calendar {
namespace internal_calendar {

class Date {
 public:
  static constexpr Date YYYYMMDD(std::int64_t const digits);
  static constexpr Date YYYYwwD(std::int64_t const digits);
  static constexpr Date YYYYDDD(std::int64_t const digits);

  static constexpr Date Calendar(int const year,
                                 int const month,
                                 int const day);
  static constexpr Date Ordinal(int const year, int const day);
  static constexpr Date Week(int const year, int const week, int const day);

  constexpr int year() const;
  constexpr int month() const;
  constexpr int day() const;

  constexpr int ordinal() const;

  constexpr Date next_day() const;

  constexpr bool operator==(Date const& other) const;
  constexpr bool operator!=(Date const& other) const;

  constexpr bool operator<(Date const& other) const;
  constexpr bool operator>(Date const& other) const;
  constexpr bool operator<=(Date const& other) const;
  constexpr bool operator>=(Date const& other) const;

 private:
  constexpr Date(int const year,
                 int const month,
                 int const day);

  int const year_;
  int const month_;
  int const day_;
};

class Time {
 public:
  static constexpr Time hhmmss_ms(int const hhmmss, int ms);

  constexpr int hour() const;
  constexpr int minute() const;
  constexpr int second() const;
  constexpr int millisecond() const;

  constexpr bool is_leap_second() const;
  // Whether |*this| is 24:00:00.
  constexpr bool is_end_of_day() const;

 private:
  constexpr Time(int const hour,
                 int const minute,
                 int const second,
                 int const millisecond);

  // Checks that this represents a valid time of day as per ISO 8601, thus
  // that the components are in the normal range, or that the object represents
  // a time in a leap second, or that it represents the end of the day.
  constexpr Time const& checked() const;

  int const hour_;
  int const minute_;
  int const second_;
  int const millisecond_;
};

class DateTime {
 public:
  static constexpr DateTime BeginningOfDay(Date const& date);

  constexpr Date const& date() const;
  constexpr Time const& time() const;

  // If |time()| is 24:00:00, returns an equivalent DateTime where midnight is
  // expressed as 00:00:00 on the next day; otherwise, returns |*this|.
  constexpr DateTime normalized_end_of_day() const;

 private:
  constexpr DateTime(Date const date, Time const time);

  // Checks that |time| does not represent a leap second unless |date| is the
  // last day of the month.
  constexpr DateTime const& checked() const;

  Date const date_;
  Time const time_;

  friend constexpr DateTime operator""_DateTime(char const* str,
                                                std::size_t size);
};

constexpr Date operator""_Date(char const* str, std::size_t size);
constexpr Time operator""_Time(char const* str, std::size_t size);
constexpr DateTime operator""_DateTime(char const* str, std::size_t size);


// Returns the duration between 2000-01-01T12:00:00 and |date_time| (of the same
// timescale), not counting any leap seconds that may have occurred in the past.
// |date_time| itself may be leap second.
// Note that this may count non-SI seconds depending on the time scale according
// to which it is interpreted.
// On a time scale with leap seconds, this is not injective: a positive leap
// second and the following second map to the same interval.
constexpr quantities::Time TimeScale(DateTime const& date_time);

}  // namespace internal_calendar

using internal_calendar::Date;
using internal_calendar::DateTime;
using internal_calendar::operator""_Date;
using internal_calendar::operator""_DateTime;
using internal_calendar::operator""_Time;
using internal_calendar::Time;
using internal_calendar::TimeScale;

}  // namespace calendar
}  // namespace astronomy
}  // namespace principia

#include "astronomy/calendar_body.hpp"
