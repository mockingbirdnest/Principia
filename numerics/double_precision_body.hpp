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

template<typename T>
void DoublePrecision<T>::WriteToMessage(
    not_null<serialization::DoublePrecision*> const message) const {
  value.WriteToMessage(message->mutable_value());
  error.WriteToMessage(message->mutable_error());
}

template<typename T>
DoublePrecision<T> DoublePrecision<T>::ReadFromMessage(
    serialization::DoublePrecision const& message) {
  value = T::ReadFromMessage(message.value());
  error = T::ReadFromMessage(message.error());
}

}  // namespace numerics
}  // namespace principia
