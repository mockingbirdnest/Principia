#pragma once

#include "geometry/interval.hpp"

#include <algorithm>

namespace principia {
namespace geometry {
namespace internal_interval {

using quantities::NaN;

template<typename T>
Difference<T> Interval<T>::measure() const {
  return max >= min ? max - min : Difference<T>{};
}

template<typename T>
T Interval<T>::midpoint() const {
  return max >= min ? min + measure() / 2 : NaN<T>;
}

template<typename T>
void Interval<T>::Include(T const& x) {
  min = std::min(min, x);
  max = std::max(max, x);
}

}  // namespace internal_interval
}  // namespace geometry
}  // namespace principia
