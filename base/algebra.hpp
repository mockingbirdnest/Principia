#pragma once

#include <concepts>
#include <type_traits>

#include "base/not_constructible.hpp"
#include "base/traits.hpp"

namespace principia {
namespace base {
namespace _algebra {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_traits;

// TODO(egg): additive_group should subsume affine, but we use it there.
template<typename G>
concept additive_group = requires(G x, G y, int n) {
  G{};
  // For some reason, unary + on a boost number returns a const number.  The
  // other operators correctly return a number (with et_off).
  { +x } -> std::convertible_to<G>;
  { -x } -> std::same_as<G>;
  { x + y } -> std::same_as<G>;
  { x - y } -> std::same_as<G>;
  { x += y } -> std::same_as<G&>;
  { x -= y } -> std::same_as<G&>;
  // An abelian group is a ℤ-module; we require the corresponding
  // operations. Note that `std::integral`, not `int`, should be used when
  // implementing these operations to avoid implicit conversions from
  // `double`.
  { n * x } -> std::same_as<G>;
  { x * n } -> std::same_as<G>;
  { x *= n } -> std::same_as<G&>;
};

// A set acted upon simply transitively by an additive group.
template<typename A>
concept affine = requires(A x, A y) {
  { x - y } -> additive_group;
  { y + (x - y) } -> std::same_as<A>;
  { y += (x - y) } -> std::same_as<A&>;
  { y - (x - y) } -> std::same_as<A>;
  { y -= (x - y) } -> std::same_as<A&>;
};


// The result type of + and - on arguments of types `Left` and `Right`.
// The operators must be arithmetic: Sum<std::string, std::string> is ill-formed.
template<affine Left, affine Right>
requires additive_group<Left> || additive_group<Right>
using Sum = decltype(std::declval<Left>() + std::declval<Right>());
template<affine Left, affine Right = Left>
using Difference = decltype(std::declval<Left>() - std::declval<Right>());

// A graded ring restricted to its homogeneous elements; multiplication can
// alter the type, and addition is only defined between homogeneous
// elements.
template<typename A>
concept homogeneous_ring = additive_group<A> && requires(A x, A y, A z) {
  // Really, multiplication should return a homogeneous_ring.
  { x * y } -> additive_group;
  { (x * y) * z } -> additive_group;
  { (x * y) * z } -> std::same_as<decltype(x * (y * z))>;
};

template<typename A>
concept ring = homogeneous_ring<A> && requires(A x, A y) {
  { x * y } -> std::same_as<A>;
  { x *= y } -> std::same_as<A&>;
} && (requires {
  { 1 } -> std::convertible_to<A>;
} || requires {
  { A::Identity() } -> std::same_as<A>;
});

// TODO(egg): field should subsume homogeneous_field, but we use it in
// homogeneous_field.

template<typename K>
concept field = ring<K> && !std::integral<K> && requires(K x, K y, K z) {
  { 1 / y } -> std::same_as<K>;
  { x / y } -> std::same_as<K>;
  { x /= y } -> std::same_as<K&>;
};

// TODO(egg): vector_space should subsume homogeneous_vector_space, but we
// use it in homogeneous_vector_space.

template<typename V, typename K>
concept vector_space = field<K> && requires(K λ, V v) {
  { λ * v } -> std::same_as<V>;
  { v * λ } -> std::same_as<V>;
  { v / λ } -> std::same_as<V>;
  { v *= λ } -> std::same_as<V&>;
  { v /= λ } -> std::same_as<V&>;
};

// A graded field restricted to its homogeneous elements; multiplication and
// division can alter the type, and addition is only defined between
// homogeneous elements.
template<typename K>
concept homogeneous_field = homogeneous_ring<K> && requires(K x, K y, K z) {
  { x / y } -> field;
  requires vector_space<K, decltype(x / y)>;
  { (1 / x) } -> homogeneous_ring;
  { 1 / (x * y) } -> std::same_as<decltype((1 / x) * (1 / y))>;
  { (x * y) / z } -> std::same_as<K>;
  { (x / y) * z } -> std::same_as<K>;
};

template<typename V, typename K>
concept homogeneous_vector_space =
    additive_group<V> && homogeneous_field<K> && requires(K λ, K μ, V v) {
      // Really, these operations should return a homogeneous_vector_space.
      requires vector_space<V, decltype(λ / μ)>;
      { λ * v } -> additive_group;
      { v * λ } -> std::same_as<decltype(λ * v)>;
      { v / λ } -> additive_group;
      { λ * v / λ } -> std::same_as<V>;
      { (λ * μ) * v } -> additive_group;
      { λ * (μ * v) } -> std::same_as<decltype((λ * μ) * v)>;
    };

template<typename V>
concept real_vector_space = vector_space<V, double>;

template<typename A, typename K>
concept affine_space = affine<A> && requires(A x, A y) {
  { x - y } -> vector_space<K>;
};

template<typename V>
concept real_affine_space = affine_space<V, double>;

template<typename T1, typename T2>
concept hilbert =
    requires(T1 const& t1, T2 const& t2) { InnerProduct(t1, t2); };

// The result type of * and / on arguments of types `Left` and `Right`.
// TODO(egg): What requirements would be appropriate for these?
template<typename Left, typename Right>
using Product = decltype(std::declval<Left>() * std::declval<Right>());
template<typename Left, typename Right>
using Quotient = decltype(std::declval<Left>() / std::declval<Right>());

template<typename Q>
using Inverse = Quotient<double, Q>;

template<typename T, int exponent>
struct ExponentiationGenerator;

// The type of iterated multiplication or division.
// TODO(egg): What requirements would be appropriate for these? The
// implementation certainly requires that x / y / z have the same type as
// x / (y * z).
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

using internal::additive_group;
using internal::affine;
using internal::affine_space;
using internal::field;
using internal::hilbert;
using internal::homogeneous_field;
using internal::homogeneous_ring;
using internal::homogeneous_vector_space;
using internal::real_affine_space;
using internal::real_vector_space;
using internal::ring;
using internal::vector_space;
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

}  // namespace _algebra

}  // namespace base
}  // namespace principia

#include "base/algebra_body.hpp"
