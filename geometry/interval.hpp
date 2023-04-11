#pragma once

#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _interval {
namespace internal {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Represents the interval [min, max]; T must be an ordered affine space.
template<typename T>
struct Interval {
  T min = T{} + Infinity<Difference<T>>;  // NOLINT(whitespace/braces)
  T max = T{} - Infinity<Difference<T>>;  // NOLINT(whitespace/braces)

  // The Lebesgue measure of this interval.
  Difference<T> measure() const;
  // The midpoint of this interval; NaN if the interval is empty (min > max).
  T midpoint() const;

  // Extends this interval so that it contains x.
  void Include(T const& x);
};

template<typename T>
std::ostream& operator<<(std::ostream& out, Interval<T> const& interval);

}  // namespace internal

using internal::Interval;

}  // namespace _interval
}  // namespace geometry
}  // namespace principia

#include "geometry/interval_body.hpp"
