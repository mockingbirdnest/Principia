#pragma once

#include "base/algebra.hpp"

namespace principia {
namespace base {
namespace _algebra {
namespace internal {

template<typename T, int exponent>
  requires(exponent > 1)
struct ExponentiationGenerator<T, exponent> : not_constructible {
  using type =
      Product<T, typename ExponentiationGenerator<T, exponent - 1>::type>;
};
template<typename T, int exponent>
  requires(exponent <= 0)
struct ExponentiationGenerator<T, exponent> : not_constructible {
  using type =
      Quotient<T, typename ExponentiationGenerator<T, -exponent + 1>::type>;
};
template<typename T>
struct ExponentiationGenerator<T, 1> : not_constructible {
  using type = T;
};

}  // namespace internal
}  // namespace _algebra
}  // namespace base
}  // namespace principia
