
#pragma once

#include "geometry/named_quantities.hpp"

namespace principia {
namespace astronomy {
namespace internal_date {

using geometry::Instant;

constexpr Instant operator""_TAI(char const* string, std::size_t size);
constexpr Instant operator""_TT(char const* string, std::size_t size);
constexpr Instant operator""_UTC(char const* string, std::size_t size);

}  // namespace internal_date

using internal_date::operator""_TAI;
using internal_date::operator""_TT;
using internal_date::operator""_UTC;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/date_body.hpp"
