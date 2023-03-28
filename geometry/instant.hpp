#pragma once

#include <ostream>
#include <string>

#include "geometry/point.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _instant {
namespace internal {

using namespace principia::geometry::_point;
using namespace principia::quantities::_quantities;

using Instant = Point<Time>;

constexpr Instant InfinitePast = Instant() - Infinity<Time>;
constexpr Instant InfiniteFuture = Instant() + Infinity<Time>;

// IEEE 754:2008 nextUp and nextDown for Instants.
// We would like to avoid the terms “up” and “down” when referring to the
// passage of time.  We avoid the term “next” in one direction because of the
// confusability with |std::nextafter|, which has different semantics, and in
// the other because of the awkwardness of the phrase “next before”.
// Defined inline for want of a way to alias functions in C++.
constexpr Instant JustAfter(Instant const& t) { return NextUp(t); }
constexpr Instant JustBefore(Instant const& t) { return NextDown(t); }

std::string DebugString(Instant const& t);
std::ostream& operator<<(std::ostream& os, Instant const& t);

}  // namespace internal

using internal::InfiniteFuture;
using internal::InfinitePast;
using internal::Instant;
using internal::JustAfter;
using internal::JustBefore;
using internal::operator<<;

}  // namespace _instant
}  // namespace geometry
}  // namespace principia
