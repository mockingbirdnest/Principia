
#pragma once

#include "geometry/hilbert.hpp"

#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {
namespace internal_hilbert {

template<typename T1, typename T2, typename U>
auto Hilbert<T1, T2, U>::InnerProduct(T1 const& t1, T2 const& t2)
    -> InnerProductType {
  return t1 * t2;
}

template<typename T1, typename T2>
auto Hilbert<T1, T2,
             std::void_t<decltype(InnerProduct(std::declval<T1>(),
                                               std::declval<T2>()))>>::
    InnerProduct(T1 const& t1, T2 const& t2) -> InnerProductType {
  // Is there a better way to avoid recursion than to put our fingers inside
  // grassmann?
  return internal_grassmann::InnerProduct(t1, t2);
}

}  // namespace internal_hilbert
}  // namespace geometry
}  // namespace principia
