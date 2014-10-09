#pragma once

#include "geometry/named_quantities.hpp"
#include "quantities/si.hpp"

// Epochs for expressing instants of Temps Atomique International (or Terrestial
// Time, which is parallel to it).
namespace principia {
namespace geometry {

// |kJ2000| represents to the standard epoch J2000.0.
// According to Resolution B1 (On the Use of Julian Dates) of the XXIIIrd IAU
// general assembly, "it is recommended that JD be specified as SI seconds in
// Terrestrial Time (TT)", see http://goo.gl/oPemRm. J2000.0 is by definition
// JD 2451545.0, i.e., noon on the first of January, 2000 (TT). One computes the
// corresponding Temps Atomique International using TT = TAI + 32.184 s. The
// corresponding UTC is obtained by subtracting 10 s of inital offset and 22
// leap seconds.
// +2000-01-01T12:00:00.000 (Terrestial Time).
// +2000-01-01T11:59:27.816 (Temps Atomique International).
// +2000-01-01T11:58:55.816 (UTC).
Instant const kJ2000 = Instant();

// Unix epoch, obtained by subtracting 30 years (including 7 leap years)
// and the UTC time of day at J2000.0 from J2000.0.
// +2000-01-01T00:00:00.000 (UTC)
Instant const kUnixEpoch = kJ2000 - ((30 * 365 + 7) * si::Day + 11 * si::Hour +
                                     58 * si::Minute + 55.816 * si::Second);

// The Julian Date JD |days|. J2000.0 is JD 2451545.0. |days| is the number of
// days since -4712-01-01-T12:00:00.000 (Terrestrial Time, Julian calendar).
Instant JulianDate(double const days);
// The Modified Julian Date MJD |days|. MJD is defined as JD - 2400000.5 days,
// so |ModifiedJulianDate(0)| is +1858-11-17-T00:00:00.000 (Terrestrial Time).
Instant ModifiedJulianDate(double const days);

}  // namespace geometry
}  // namespace principia

#include "geometry/epoch_body.hpp"
