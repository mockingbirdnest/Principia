#pragma once

#include "base/macros.hpp"

namespace principia {
namespace integrators {

template<typename Scalar>
inline DoublePrecision<Scalar>::DoublePrecision(Scalar const& value)
    : value(value) {}

template<typename Scalar>
FORCE_INLINE DoublePrecision<Scalar>::Increment(Scalar const& increment) {
  // The naming conventions follow Higham, Accuracy and Stability of Numerical
  // Algorithms, Algorithm 4.2.
  Scalar const temp = value;
  Scalar const y = increment + error;
  value = temp + y;
  error = (temp - value) + y;
}

}  // namespace integrators
}  // namespace principia
