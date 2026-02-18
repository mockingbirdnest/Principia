#pragma once

#include "geometry/hilbert.hpp"

#include "numerics/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _hilbert {
namespace internal {

using namespace principia::numerics::_elementary_functions;

template<homogeneous_field T, homogeneous_field U>
  requires std::totally_ordered<Product<T, U>>
constexpr Product<T, U> InnerProduct(T const& left, U const& right) {
  return left * right;
}

template<typename T>
  requires (requires(T x) { x.Norm²(); })
constexpr decltype(std::declval<T>().Norm²()) Norm²(T const& x) {
  return x.Norm²();
}

template<typename T>
  requires (requires(T x) { x.Norm(); })
constexpr decltype(std::declval<T>().Norm()) Norm(T const& x) {
  return x.Norm();
}

template<homogeneous_field T>
  requires std::totally_ordered<T>
constexpr T Norm(T const& x) {
  return Abs(x);
}

}  // namespace internal
}  // namespace _hilbert
}  // namespace geometry
}  // namespace principia
