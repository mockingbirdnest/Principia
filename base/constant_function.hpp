#pragma once

namespace principia {
namespace base {

// A function object that can be called with arbitrary arguments, and always
// returns the same |value|.
template<typename T>
struct ConstantFunction {
  template<typename... Args>
  constexpr T operator()(Args&&...) {
    return value;
  }

  T value;
};

template <typename T>
constexpr ConstantFunction<T> Identically(T value) {
  return {value};
}

}  // namespace base
}  // namespace principia
