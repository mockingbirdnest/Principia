
#pragma once

#include "geometry/named_quantities.hpp"
#include "quantities/si.hpp"

// |geometry::Instant| represents instants of Terrestrial Time (TT).  The
// utilities in this file provide its standard epoch and two ways of specifying
// TT dates.

namespace principia {

using geometry::Instant;
using quantities::si::Second;

namespace astronomy {

// |J2000| represents to the standard epoch J2000.0.
// According to Resolution B1 (On the Use of Julian Dates) of the XXIIIrd IAU
// general assembly, "it is recommended that JD be specified as SI seconds in
// Terrestrial Time (TT)", see http://goo.gl/oPemRm. J2000.0 is by definition
// JD 2451545.0, i.e., noon on the first of January, 2000 (TT).
// "2000-01-01T12:00:00"_TT
// "2000-01-01T11:59:27,816"_TAI
// "2000-01-01T11:58:55,816"_UTC
constexpr Instant J2000;

constexpr Instant InfinitePast =
    J2000 - std::numeric_limits<double>::infinity() * Second;
constexpr Instant InfiniteFuture =
    J2000 + std::numeric_limits<double>::infinity() * Second;

// The Julian Date JD |days|. J2000.0 is JD 2451545.0. |days| is the number of
// days since -4712-01-01T12:00:00,000 (Terrestrial Time, Julian calendar).
constexpr Instant JulianDate(double const days);
// The Modified Julian Date MJD |days|. MJD is defined as JD - 2400000.5 days,
// so |ModifiedJulianDate(0)| is "1858-11-17T00:00:00"_TT.
constexpr Instant ModifiedJulianDate(double const days);

}  // namespace astronomy
}  // namespace principia

#include "astronomy/epoch_body.hpp"
