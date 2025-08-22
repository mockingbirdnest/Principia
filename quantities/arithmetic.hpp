#pragma once
#pragma once

#include "base/not_constructible.hpp"

// This file defines alias templates that are derived from the operators +, -,
// *, and /.  These are particularly useful with quantities and types built on
// top of `Quantity`, but they are not specific to `Quantity`, and have no
// dependency on anything else in the project `quantities`.

namespace principia {
namespace quantities {
namespace _arithmetic {
namespace internal {

using namespace principia::base::_not_constructible;

// The result type of +, -, * and / on arguments of types `Left` and `Right`.
template<typename Left, typename Right>
using Sum = decltype(std::declval<Left>() + std::declval<Right>());
template<typename Left, typename Right = Left>
using Difference = decltype(std::declval<Left>() - std::declval<Right>());
template<typename Left, typename Right>
using Product = decltype(std::declval<Left>() * std::declval<Right>());
template<typename Left, typename Right>
using Quotient = decltype(std::declval<Left>() / std::declval<Right>());

template<typename Q>
using Inverse = Quotient<double, Q>;

template<typename T, int exponent>
struct ExponentiationGenerator;

// The type of iterated multiplication or iterated
template<typename T, int exponent>
using Exponentiation = typename ExponentiationGenerator<T, exponent>::type;
template<typename Q>
using Square = Exponentiation<Q, 2>;
template<typename Q>
using Cube = Exponentiation<Q, 3>;

// The result type of the N-th derivative of a `Value`-valued function with
// respect to its `Argument`-valued argument.
template<typename Value, typename Argument, int order = 1>
using Derivative = typename std::conditional_t<
    order == 0,
    Value,
    Quotient<Difference<Value>, Exponentiation<Difference<Argument>, order>>>;

// The result type of the primitive of a `Value`-valued function with respect to
// its `Argument`-valued argument.  The primitive of an affine-valued function
// does not make much sense, but it must compile, hence the Difference.
template<typename Value, typename Argument>
using Primitive = Product<Difference<Value>, Difference<Argument>>;

}  // namespace internal

using internal::Cube;
using internal::Derivative;
using internal::Difference;
using internal::Exponentiation;
using internal::Inverse;
using internal::Primitive;
using internal::Product;
using internal::Quotient;
using internal::Square;
using internal::Sum;

}  // namespace _arithmetic
}  // namespace quantities
}  // namespace principia

#include "quantities/arithmetic_body.hpp"
