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
constexpr auto Norm²(T const& x) {
  return x.Norm²();
}

template<homogeneous_field T>
constexpr auto Norm²(T const& x) {
  return x * x;
}

template<typename T>
constexpr auto Norm(T const& x) {
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
