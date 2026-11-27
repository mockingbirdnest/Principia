#pragma once

#include <concepts>
#include <type_traits>

#include "base/algebra.hpp"

namespace principia {
namespace geometry {
namespace _hilbert {
namespace internal {

using namespace principia::base::_algebra;

template<homogeneous_field T, homogeneous_field U>
  requires std::totally_ordered<Product<T, U>>
constexpr Product<T, U> InnerProduct(T const& left, U const& right);

template<typename T, typename U = T>
using InnerProductType =
    decltype(InnerProduct(std::declval<T>(), std::declval<U>()));

template<typename T>
  requires(requires(T x) { x.Norm²(); })
constexpr auto Norm²(T const& x);

template<homogeneous_field T>
  requires std::totally_ordered<T>
constexpr Square<T> Norm²(T const& x);

template<typename T>
using Norm²Type = decltype(Norm²(std::declval<T>()));

template<typename T>
  requires(requires(T x) { x.Norm(); })
constexpr auto Norm(T const& x);

template<homogeneous_field T>
  requires std::totally_ordered<T>
constexpr T Norm(T const& x);

template<typename T>
using NormType = decltype(Norm(std::declval<T>()));

template<typename T>
using NormalizedType = Quotient<T, NormType<T>>;

template<typename T>
  requires(requires(T x) { Norm(x); })
constexpr NormalizedType<T> Normalize(T const& x);

template<typename V>
concept hilbert = requires(V u, V v) {
  { InnerProduct(u, v) } -> homogeneous_field;
  { Norm(u) } -> homogeneous_field;
  { Norm²(u) } -> std::same_as<Square<NormType<V>>>;
  // These types need not be the same, e.g., in ℂ, the InnerProductType is ℂ but
  // the Norm²Type is ℝ.  Norm²Type should be a subfield of InnerProductType.
  requires std::convertible_to<Norm²Type<V>, InnerProductType<V>>;
  requires std::totally_ordered<Norm²Type<V>>;
  requires std::totally_ordered<NormType<V>>;
  requires homogeneous_vector_space<V, NormType<V>>;
};

template<typename A>
concept euclidean = affine<A> && hilbert<Difference<A>>;

template<typename A, typename K>
concept euclidean_space = affine_space<A, K> && hilbert<Difference<A>>;

template<typename A, typename K>
concept homogeneous_euclidean_space =
    homogeneous_affine_space<A, K> && hilbert<Difference<A>>;

}  // namespace internal

using internal::euclidean;
using internal::euclidean_space;
using internal::hilbert;
using internal::homogeneous_euclidean_space;
using internal::InnerProduct;
using internal::InnerProductType;
using internal::Norm;
using internal::Normalize;
using internal::NormalizedType;
using internal::NormType;
using internal::Norm²;
using internal::Norm²Type;

}  // namespace _hilbert
}  // namespace geometry
}  // namespace principia

#include "geometry/hilbert_body.hpp"
