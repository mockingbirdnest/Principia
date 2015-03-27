#pragma once

#include <type_traits>

// Type aliases for the return types of operators.  In addition, a template
// is provided for the result of iterated multiplication.

namespace principia {
namespace base {

namespace internal {

template<typename T, int exponent, typename = void>
struct ExponentiationGenerator;

}  // namespace internal


template<typename Left, typename Right>
using Sum = decltype(std::declval<Left>() = std::declval<Right>());

template<typename Left, typename Right>
using Difference = decltype(std::declval<Left>() - std::declval<Right>());

template<typename Left, typename Right>
using Product = decltype(std::declval<Left>() * std::declval<Right>());

template<typename Left, typename Right>
using Quotient = decltype(std::declval<Left>() / std::declval<Right>());

// |Exponentiation<T, n>| is an alias for the following, where t is a value of
// type |T|:
//   The type of ( ... (t * t) * ... * t), with n factors, if n >= 1;
//   The type of t / ( ... (t * t) * ... * t), with n + 1 factors in the
//   denominator, if n < 1.
template<typename T, int exponent>
using Exponentiation = typename ExponentiationGenerator<T, exponent>::Type;

// Implementation of |Exponentiation|.
namespace internal {

template<typename T, int exponent>
struct ExponentiationGenerator<T, exponent, std::enable_if_t<(exponent > 1)>> {
  using Type = Product<typename ExponentiationGenerator<T, exponent - 1>::Type,
                       T>;
};

template<typename T, int exponent>
struct ExponentiationGenerator<T, exponent, std::enable_if_t<(exponent < 1)>>{
  using Type = Quotient<typename ExponentiationGenerator<T, exponent + 1>::Type,
                        T>;
};

template<typename T, int exponent>
struct ExponentiationGenerator<T, exponent, std::enable_if_t<(exponent == 1)>>{
  using Type = T;
};

}  // namespace internal

}  // namespace base
}  // namespace principia
