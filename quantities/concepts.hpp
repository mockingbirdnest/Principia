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

template<typename G>
concept additive_group = requires(G x, G y) {
  G{};
  { +x } -> std::same_as<G>;
  { -x } -> std::same_as<G>;
  { x + y } -> std::same_as<G>;
  { x - y } -> std::same_as<G>;
  { x += y } -> std::same_as<G&>;
  { x -= y } -> std::same_as<G&>;
};

template<typename V, typename K>
concept vector_space = additive_group<V> && requires(K λ, V v) {
  { λ * v } -> std::same_as<V>;
  { v * λ } -> std::same_as<V>;
  { v / λ } -> std::same_as<V>;
  { v *= λ } -> std::same_as<V&>;
  { v /= λ } -> std::same_as<V&>;
};

template<typename V>
concept real_vector_space = vector_space<V, double>;

template<typename A, typename K>
concept affine_space = requires(A x, A y) {
  { x - y } -> vector_space<K>;
  { y + (x - y) } -> std::same_as<A>;
  { y += (x - y) } -> std::same_as<A&>;
  { y - (x - y) } -> std::same_as<A>;
  { y -= (x - y) } -> std::same_as<A&>;
};

template<typename V>
concept real_affine_space = affine_space<V, double>;

// std::integral || std::floating_point rather than
// std::convertible_to<double, T> because
// the former introduces ambiguities on Sign * Vector.
template<typename T>
concept convertible_to_quantity =
    std::integral<T> || std::floating_point<T> ||
    is_instance_of_v<Quantity, std::remove_cvref_t<T>>;

}  // namespace internal

using internal::additive_group;
using internal::affine_space;
using internal::convertible_to_quantity;
using internal::real_affine_space;
using internal::real_vector_space;
using internal::vector_space;

}  // namespace _concepts
}  // namespace quantities
}  // namespace principia
