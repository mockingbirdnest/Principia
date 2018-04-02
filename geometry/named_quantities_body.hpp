
#pragma once

#include "geometry/named_quantities.hpp"

namespace principia {
namespace geometry {

template<typename T, typename Void>
typename Normed<T, Void>::NormType Normed<T, Void>::Norm(T const& vector) {
  return quantities::Abs(vector);
}

template<typename T>
typename Normed<T, std::void_t<decltype(std::declval<T>().Norm())>>::NormType
Normed<T, std::void_t<decltype(std::declval<T>().Norm())>>::Norm(
    T const& vector) {
  return vector.Norm();
}

}  // namespace geometry
}  // namespace principia
