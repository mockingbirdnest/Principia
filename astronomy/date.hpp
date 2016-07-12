
#pragma once

#include "geometry/named_quantities.hpp"

namespace principia {
namespace astronomy {
namespace internal_date {

using geometry::Instant;

// NOTE(egg): We cannot use literal operator templates for strings, so if an
// invalid date is given and the result does not need to be constexpr the
// literal will fail at runtime; if proposal N3599 ever gets anywhere we'll be
// able to solve this.
// See http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3599.html,
// http://wg21.cmeerw.net/ewg/issue66.
// FWIW it seems that clang supports this proposal with
// -Wno-gnu-string-literal-operator-template.

#if (PRINCIPIA_COMPILER_CLANG || PRINCIPIA_COMPILER_CLANG_CL) && WE_LIKE_N3599
template<typename C, C... str>
constexpr Instant operator""_TAI();
template<typename C, C... str>
constexpr Instant operator""_TT();
template<typename C, C... str>
constexpr Instant operator""_UTC();
template<typename C, C... str>
constexpr Instant operator""_UT1();
#else
constexpr Instant operator""_TAI(char const* str, std::size_t size);
constexpr Instant operator""_TT(char const* str, std::size_t size);
constexpr Instant operator""_UTC(char const* str, std::size_t size);
constexpr Instant operator""_UT1(char const* str, std::size_t size);
#endif

}  // namespace internal_date

using internal_date::operator""_TAI;
using internal_date::operator""_TT;
using internal_date::operator""_UTC;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/date_body.hpp"
