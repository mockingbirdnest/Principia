
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

constexpr Angle EarthRotationAngle(Instant const tt);

constexpr Instant operator""_TAI(char const* str, std::size_t size);
constexpr Instant operator""_TT(char const* str, std::size_t size);
constexpr Instant operator""_UTC(char const* str, std::size_t size);
constexpr Instant operator""_UT1(char const* str, std::size_t size);

Instant ParseTAI(std::string const& s);
Instant ParseTT(std::string const& s);
Instant ParseUTC(std::string const& s);
Instant ParseUT1(std::string const& s);

}  // namespace internal_time_scales

using internal_time_scales::EarthRotationAngle;
using internal_time_scales::operator""_TAI;
using internal_time_scales::operator""_TT;
using internal_time_scales::operator""_UTC;
using internal_time_scales::operator""_UT1;
using internal_time_scales::ParseTAI;
using internal_time_scales::ParseTT;
using internal_time_scales::ParseUTC;
using internal_time_scales::ParseUT1;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/time_scales_body.hpp"
