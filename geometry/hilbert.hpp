#pragma once

#include <concepts>
#include <type_traits>

#include "base/algebra.hpp"

namespace principia {
namespace geometry {
namespace _hilbert {
namespace internal {

using namespace principia::base::_algebra;

template<typename T>
struct Norm²TypeGenerator {
  using type = decltype(std::declval<T>().Norm²());
};

template<homogeneous_field T>
  requires std::totally_ordered<T>
struct Norm²TypeGenerator<T> {
  using type = Square<T>;
};

template<typename T>
using Norm²Type = typename Norm²TypeGenerator<T>::type;

template<homogeneous_field T, homogeneous_field U>
  requires std::totally_ordered<Product<T, U>>
constexpr Product<T, U> InnerProduct(T const& left, U const& right);

template<typename T>
constexpr Norm²Type<T> Norm²(T const& x);

template<homogeneous_field T>
  requires std::totally_ordered<T>
constexpr Norm²Type<T> Norm²(T const& x);

template<typename T>
constexpr auto Norm(T const& x);

template<homogeneous_field T>
  requires std::totally_ordered<T>
constexpr T Norm(T const& x);

template<typename T, typename U = T>
using InnerProductType =
    decltype(InnerProduct(std::declval<T>(), std::declval<U>()));

// NOTE(egg): This is not defined in terms of Norm² because that one returns
// auto, and MSVC doesn’t like it when we declare a Norm²() that returns a
// Norm²Type computed from (a different) auto-valued Norm²().

template<typename T>
using NormType = decltype(Norm(std::declval<T>()));

template<typename T>
using NormalizedType = Quotient<T, NormType<T>>;

template<typename V>
concept hilbert = requires(V u, V v) {
  { Norm(u) } -> homogeneous_field;
  { InnerProduct(u, v) } -> std::same_as<Square<NormType<V>>>;
  { Norm²(u) } -> std::same_as<Square<NormType<V>>>;
  requires std::totally_ordered<Norm²Type<V>>;
  requires std::totally_ordered<NormType<V>>;
  requires homogeneous_vector_space<V, NormType<V>>;
};

}  // namespace internal

using internal::hilbert;
using internal::InnerProduct;
using internal::InnerProductType;
using internal::Norm;
using internal::NormalizedType;
using internal::NormType;
using internal::Norm²;
using internal::Norm²Type;

}  // namespace _hilbert
}  // namespace geometry
}  // namespace principia

#include "geometry/hilbert_body.hpp"
