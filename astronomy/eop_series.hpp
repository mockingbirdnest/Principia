
#pragma once

#include <array>

#include "astronomy/date.hpp"

namespace principia {
namespace astronomy {

// TODO(egg): we define |Date| etc. in an internal namespace, but we end up
// needing it here too; reusing internal namespaces is evil, perhaps Date should
// be exported.
namespace internal_date {

// An entry in the Experimental EOP C02 time series; represents UT1 - TAI at the
// given (UT1) |mjd|.
struct MJDToUT1MinusTAI {
  constexpr MJDToUT1MinusTAI(double const ut1_mjd,
                             quantities::Time const ut1_minus_tai);

  double const ut1_mjd;
  quantities::Time const ut1_minus_tai;
};

constexpr MJDToUT1MinusTAI::MJDToUT1MinusTAI(
    double const ut1_mjd,
    quantities::Time const ut1_minus_tai)
    : ut1_mjd(ut1_mjd),
      ut1_minus_tai(ut1_minus_tai) {}

// An entry in the EOP (IERS) 08 C04 time series; represents UT1 - UTC at
// 00:00:00 on the given |utc_date|.  The date is given as an integer of the
// form YYYY'MM'DD, which is then interpreted on demand, in order to limit the
// number of constexpr steps.
struct UTCToUT1MinusUTC {

  constexpr UTCToUT1MinusUTC(int const utc_date,
                             quantities::Time const& ut1_minus_utc);

  constexpr DateTime utc() const {
    return DateTime::BeginningOfDay(Date::YYYYMMDD(utc_date));
  }

  constexpr quantities::Time ut1() const {
    return TimeScale(utc()) + ut1_minus_utc;
  }

  constexpr quantities::Time ut1_minus_tai() const {
    return utc().date().year() >= 1972
               ? ut1_minus_utc + ModernUTCMinusTAI(utc().date())
               : ut1_minus_utc - TAIMinusStretchyUTC(utc());
  }

  int const utc_date;
  quantities::Time const ut1_minus_utc;
};

constexpr UTCToUT1MinusUTC::UTCToUT1MinusUTC(
    int const utc_date,
    quantities::Time const& ut1_minus_utc)
    : utc_date(utc_date),
      ut1_minus_utc(ut1_minus_utc) {}

#include "astronomy/experimental_eop_c02.generated.h"
#include "astronomy/eop_c04.generated.h"

// Returns the last entry in [begin, begin + size[ whose UT1 is less than or
// equal to the given |ut1|.  The range [begin, begin + size[ must be sorted by
// UT1.
// We have to use |begin| and |size| rather than |begin| and |end| because
// otherwise the MSVC complains about undefinedness of |end - begin| (even
// though Intellisense is fine with it).
constexpr MJDToUT1MinusTAI const* LookupUT1(quantities::Time const& ut1,
                                            MJDToUT1MinusTAI const* begin,
                                            std::ptrdiff_t const size) {
  return CHECKING(size > 0,
                  size == 1
                      ? CHECKING(begin->ut1_mjd <= mjd(ut1), begin)
                      : (begin + size / 2)->ut1_mjd <= mjd(ut1)
                            ? LookupUT1(ut1, begin + size / 2, size - size / 2)
                            : LookupUT1(ut1, begin, size / 2));
}

// Returns the last entry in [begin, begin + size[ whose UT1 = UTC + (UT1 - UTC)
// is less than or equal to the given |ut1|.  The range [begin, begin + size[
// must be sorted by UT1.
constexpr UTCToUT1MinusUTC const* LookupUT1(quantities::Time const& ut1,
                                            UTCToUT1MinusUTC const* begin,
                                            std::ptrdiff_t const size) {
  return CHECKING(
      size > 0,
      size == 1
          ? CHECKING(begin->ut1() <= ut1, begin)
          : (begin + size / 2)->ut1() <= ut1
                ? LookupUT1(ut1, begin + size / 2, size - size / 2)
                : LookupUT1(ut1, begin, size / 2));
}

constexpr auto lolcat = eop_c04[142];
static_assert(lolcat.utc_date==1962'05'23, "");

// Linear interpolation on the UT1 range [low->ut1(), (low + 1)->ut1()].  Note
// that we cannot use |Barycentre|, because it uses non-constexpr |std::vector|.
constexpr Instant InterpolatedEOPC04(UTCToUT1MinusUTC const* low,
                                     quantities::Time const& ut1) {
  // TODO(egg): figure out whether using the divided difference of the
  // |p->ut1_minus_tai()|s leads to less catastrophic cancellation than using
  // the divided difference of the |DateTimeAsUTC(p->utc())|s.
  return FromTAI(ut1 - (low->ut1_minus_tai() +
                        (ut1 - low->ut1()) * ((low + 1)->ut1_minus_tai() -
                                              low->ut1_minus_tai()) /
                            ((low + 1)->ut1() - low->ut1())));
}

constexpr Instant InterpolatedExperimentalEOPC02(MJDToUT1MinusTAI const* low,
                                                 quantities::Time const& ut1) {
  return FromTAI(ut1 - (low->ut1_minus_tai +
                        (mjd(ut1) - low->ut1_mjd) *
                            ((low + 1)->ut1_minus_tai - low->ut1_minus_tai) /
                            ((low + 1)->ut1_mjd - low->ut1_mjd)));
}

// Linear interpolation in the segment between the UT1s |low->ut1_mjd| and
// |eop_c04[0].ut1()|, used to get continuity when switching between the series.
constexpr Instant ExperimentalEOPC02ToEOPC04(MJDToUT1MinusTAI const* low,
                                             quantities::Time const& ut1) {
  return FromTAI(ut1 - (low->ut1_minus_tai +
                        (mjd(ut1) - low->ut1_mjd) *
                            (eop_c04[0].ut1_minus_tai() - low->ut1_minus_tai) /
                            ((mjd(eop_c04[0].ut1()) - low->ut1_mjd))));
}

constexpr Instant DateTimeAsUT1(DateTime const& ut1) {
  return TimeScale(ut1) < eop_c04[0].ut1()
             ? ((LookupUT1(TimeScale(ut1),
                           &experimental_eop_c02[0],
                           experimental_eop_c02.size()) +
                 1)->ut1_mjd > mjd(eop_c04[0].ut1())
                    ? ExperimentalEOPC02ToEOPC04(
                          LookupUT1(TimeScale(ut1),
                                    &experimental_eop_c02[0],
                                    experimental_eop_c02.size()),
                          TimeScale(ut1))
                    : InterpolatedExperimentalEOPC02(
                          LookupUT1(TimeScale(ut1),
                                    &experimental_eop_c02[0],
                                    experimental_eop_c02.size()),
                          TimeScale(ut1)))
             : InterpolatedEOPC04(
                   LookupUT1(TimeScale(ut1), &eop_c04[0], eop_c04.size()),
                   TimeScale(ut1));
}

constexpr auto this_is_an_identifier =
    DateTimeAsUT1("1964-03-31T23:59:59"_DateTime) -
    "1964-03-31T23:59:59.160"_UTC;

constexpr auto const* foo =
    LookupUT1(-1'000 * Second, &eop_c04[0], eop_c04.size());

constexpr auto const* bar = LookupUT1(-1'000 * Second,
                                      &experimental_eop_c02[0],
                                      experimental_eop_c02.size());

}  // namespace internal_date
}  // namespace astronomy
}  // namespace principia
