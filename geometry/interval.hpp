#pragma once

#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_interval {

using quantities::Difference;
using quantities::Infinity;

// Represents the interval [min, max]; T must be an ordered affine space.
template<typename T>
struct Interval {
  T min = T{} + Infinity<Difference<T>>();
  T max = T{} - Infinity<Difference<T>>();

  // The Lebesgue measure of this interval.
  Difference<T> measure() const;
  // The midpoint of this interval; NaN if the interval is empty (min > max).
  T midpoint() const;

  // Extends this interval so that it contains x.
  void Include(T const& x);
};

}  // namespace internal_interval

using internal_interval::Interval;

}  // namespace geometry
}  // namespace principia

#include "geometry/interval_body.hpp"
