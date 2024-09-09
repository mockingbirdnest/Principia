#pragma once

#include <concepts>
#include <type_traits>

#include "base/traits.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace _concepts {
namespace internal {

using namespace base::_traits;
using namespace quantities::_quantities;

// TODO(egg): additive_group should subsume affine, but we use it there.
// We use `convertible_to` here because we want this concept to work with
// Boost multiprecision types which heavily use implicit conversions.
template<typename G>
concept additive_group = requires(G x, G y, int n) {
  G{};
  { +x } -> std::convertible_to<G>;
  { -x } -> std::convertible_to<G>;
  { x + y } -> std::convertible_to<G>;
  { x - y } -> std::convertible_to<G>;
  { x += y } -> std::convertible_to<G&>;
  { x -= y } -> std::convertible_to<G&>;
  // An abelian group is a ℤ-module; we require the corresponding operations.
  { n * x } -> std::convertible_to<G>;
  { x * n } -> std::convertible_to<G>;
  { x *= n } -> std::convertible_to<G&>;
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

// A graded ring restricted to its homogeneous elements; multiplication can
// alter the type, and addition is only defined between homogeneous elements.
template<typename A>
concept homogeneous_ring =
    additive_group<A> && requires(A x, A y, A z) {
      // Really, multiplication should return a homogeneous_ring.
      { x * y } -> additive_group;
      { (x * y) * z } -> additive_group;
      { (x * y) * z } -> std::same_as<decltype(x * (y * z))>;
};

template<typename A>
concept ring = homogeneous_ring<A> && requires(A x, A y) {
  { x * y } -> std::same_as<A>;
  { x *= y } -> std::same_as<A&>;
};

// TODO(egg): field should subsume homogeneous_field, but we use it in
// homogeneous_field.

template<typename K>
concept field = ring<K> && !std::integral<K> && requires(K x, K y, K z) {
  { 1 } -> std::convertible_to<K>;
  { 1 / y } -> std::same_as<K>;
  { x / y } -> std::same_as<K>;
  { x /= y } -> std::same_as<K&>;
};

// TODO(egg): vector_space should subsume homogeneous_vector_space, but we use
// it in homogeneous_vector_space.

template<typename V, typename K>
concept vector_space = field<K> && requires(K λ, V v) {
      { λ * v } -> std::same_as<V>;
      { v * λ } -> std::same_as<V>;
      { v / λ } -> std::same_as<V>;
      { v *= λ } -> std::same_as<V&>;
      { v /= λ } -> std::same_as<V&>;
    };

// A graded field restricted to its homogeneous elements; multiplication and
// division can alter the type, and addition is only defined between homogeneous
// elements.
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

template<typename T>
concept quantity = instance<T, Quantity> || std::same_as<T, double>;

// std::integral || std::floating_point rather than
// std::convertible_to<double, T> because
// the latter introduces ambiguities on Sign * Vector.
template<typename T>
concept convertible_to_quantity =
    quantity<std::remove_cvref_t<T>> ||
    std::integral<std::remove_cvref_t<T>> ||
    std::floating_point<std::remove_cvref_t<T>>;

}  // namespace internal

using internal::additive_group;
using internal::affine;
using internal::affine_space;
using internal::convertible_to_quantity;
using internal::field;
using internal::homogeneous_field;
using internal::homogeneous_ring;
using internal::homogeneous_vector_space;
using internal::quantity;
using internal::real_affine_space;
using internal::real_vector_space;
using internal::ring;
using internal::vector_space;

}  // namespace _concepts
}  // namespace quantities
}  // namespace principia
