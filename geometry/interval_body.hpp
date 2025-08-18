#pragma once

#include "geometry/interval.hpp"

#include <algorithm>

#include "base/concepts.hpp"
#include "glog/logging.h"

namespace principia {
namespace geometry {
namespace _interval {
namespace internal {

using namespace principia::base::_concepts;

template<typename T>
Difference<T> Interval<T>::measure() const {
  return max >= min ? max - min : Difference<T>{};
}

template<typename T>
bool Interval<T>::empty() const {
  return max <= min;
}

template<typename T>
T Interval<T>::midpoint() const {
  if constexpr (boost_cpp_number<T>) {
    DCHECK_GE(max, min);
    return min + measure() / 2;
  } else {
    return max >= min ? min + measure() / 2 : min + NaN<Difference<T>>;
  }
}

template<typename T>
void Interval<T>::Include(T const& x) {
  min = std::min(min, x);
  max = std::max(max, x);
}

template<typename T>
std::ostream& operator<<(std::ostream& out, Interval<T> const& interval) {
  return out << interval.midpoint() << " ± " << interval.measure() / 2;
}

}  // namespace internal
}  // namespace _interval
}  // namespace geometry
}  // namespace principia
