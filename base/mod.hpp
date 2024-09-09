#pragma once

namespace principia {
namespace base {
namespace _mod {
namespace internal {

// The result of `mod` has the same sign as `divisor` (same convention as Ada,
// ISO Prolog, Haskell, etc.), whereas the result of `dividend % divisor` has
// the same sign as `dividend` (like `rem` in the aforementioned languages).
constexpr int mod(int const dividend, int const divisor) {
  return ((dividend % divisor) + divisor) % divisor;
}

// Similar to Mathematica's `Mod`: the result is congruent to `dividend` modulo
// `divisor`, and lies in [offset, divisor + offset[ if divisor > 0, and in
// ]divisor + offset, offset] otherwise.
constexpr int mod(int const dividend, int const divisor, int const offset) {
  return mod(dividend - offset, divisor) + offset;
}

}  // namespace internal

using internal::mod;

}  // namespace _mod
}  // namespace base
}  // namespace principia
