
#pragma once

#include <string>

#include "geometry/named_quantities.hpp"

namespace principia {
namespace astronomy {
namespace internal_time_scales {

using geometry::Instant;
using quantities::Angle;

// NOTE(egg): We cannot use literal operator templates for strings, so if an
// invalid date is given and the result does not need to be constexpr the
// literal will fail at runtime; if proposal N3599 ever gets anywhere we'll be
// able to solve this.
// See http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3599.html,
// http://wg21.cmeerw.net/ewg/issue66.
// FWIW it seems that clang supports this proposal with
// -Wno-gnu-string-literal-operator-template.

constexpr Angle EarthRotationAngle(Instant tt);

// Astronomical time scales:
// — Temps Atomique International;
// — Temps Terrestre;
// — Coordinated Universal Time;
// — Universal Time (UT1).
// |Instant| represents TT; note that the conversion from TAI, UTC, and UT1
// converts to TT(TAI), not TT(BIPM); this may cause issues where long-term
// frequency stability and high frequency accuracy are needed, see
// https://www.bipm.org/en/bipm-services/timescales/time-ftp/ttbipm.html.

constexpr Instant operator""_TAI(char const* str, std::size_t size);
constexpr Instant operator""_TT(char const* str, std::size_t size);
constexpr Instant operator""_UTC(char const* str, std::size_t size);
constexpr Instant operator""_UT1(char const* str, std::size_t size);

Instant ParseTAI(std::string const& s);
Instant ParseTT(std::string const& s);
Instant ParseUTC(std::string const& s);
Instant ParseUT1(std::string const& s);

// GNSS time scales.
// As documented, e.g., in p. 32 of the RINEX specification, v. 3.03, update 1,
// ftp://igs.org/pub/data/format/rinex303_update1.pdf,
// apart from the small errors in the realizations of the different time
// systems,
//   TAI - 19 s = GPS Time
//              = Galileo System Time (GST)
//              = 準天頂衛星 (Quasi-Zenith Satellite) Time (QZST)
//              = IRNSS Network Time (IRNWT),
//   TAI - 33 s = 北斗 (BeiDou) Time (北斗时, BDT),
//   UTC        = ГЛОНАСС Time.
// The errors in realization are on the order of nanoseconds (nominally 1 μs for
// GPS, 50 ns for Galileo, 50 ns for 北斗).
// Note that by not distinguishing TT(TAI) from TT(BIPM), we disregard much
// larger errors between TAI and TT; for instance,
//   TT(BIPM17) = TAI + 32.184 s + 27661.0 ns
// at the end of 2017.
// Even ignoring the 27.6 μs offset, TT(BIPMyy) varies with respect to TAI on
// the order of nanoseconds over a few weeks.
// Since we disregard small errors in realization, we identify ГЛОНАСС time with
// UTC, and the GNSS time scales other than 北斗 and ГЛОНАСС with GPS time.
// We do not support Julian dates in GNSS time scales.

constexpr Instant operator""_GPS(char const* str, std::size_t size);
constexpr Instant operator""_北斗(char const* str, std::size_t size);

Instant ParseGPSTime(std::string const& s);
Instant Parse北斗Time(std::string const& s);

}  // namespace internal_time_scales

using internal_time_scales::EarthRotationAngle;
using internal_time_scales::Parse北斗Time;
using internal_time_scales::ParseGPSTime;
using internal_time_scales::ParseTAI;
using internal_time_scales::ParseTT;
using internal_time_scales::ParseUT1;
using internal_time_scales::ParseUTC;
using internal_time_scales::operator""_北斗;
using internal_time_scales::operator""_GPS;
using internal_time_scales::operator""_TAI;
using internal_time_scales::operator""_TT;
using internal_time_scales::operator""_UT1;
using internal_time_scales::operator""_UTC;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/time_scales_body.hpp"
