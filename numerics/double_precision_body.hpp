#pragma once

#include "numerics/double_precision.hpp"

#include "base/macros.hpp"

namespace principia {
namespace numerics {

template<typename T>
inline DoublePrecision<T>::DoublePrecision(T const& value)
    : value(value),
      error() {}

template<typename T>
FORCE_INLINE void DoublePrecision<T>::Increment(
    Difference<T> const& increment) {
  // The naming conventions follow Higham, Accuracy and Stability of Numerical
  // Algorithms, Algorithm 4.2.
  T const temp = value;
  Difference<T> const y = increment + error;
  value = temp + y;
  error = (temp - value) + y;
}

}  // namespace numerics
}  // namespace principia
