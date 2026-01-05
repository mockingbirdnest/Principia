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

template<affine Value, affine Argument>
struct DerivativeGenerator<Value, Argument, 0> : not_constructible {
  using type = Value;
};

template<affine Value, affine Argument, int order>
  requires(order > 0) && affine_module<Value, Difference<Argument>>
struct DerivativeGenerator<Value, Argument, order> : not_constructible {
  using type = Difference<Value>;
};

template<affine Value, affine Argument, int order>
  requires(order > 0) && !affine_module<Value, Difference<Argument>> &&
          homogeneous_affine_space<Value, Difference<Argument>>
struct DerivativeGenerator<Value, Argument, order> : not_constructible {
  using type =
      Quotient<Difference<Value>, Exponentiation<Difference<Argument>, order>>;
};

}  // namespace internal
}  // namespace _algebra
}  // namespace base
}  // namespace principia
