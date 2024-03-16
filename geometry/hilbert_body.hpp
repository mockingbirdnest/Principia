#pragma once

#include "geometry/hilbert.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _hilbert {
namespace internal {

using namespace principia::quantities::_elementary_functions;

template<typename T1, typename T2>
  requires quantity<T1> && quantity<T2> && (!std::is_same_v<T1, T2>)
auto Hilbert<T1, T2>::InnerProduct(T1 const& t1, T2 const& t2)
    -> InnerProductType {
  return t1 * t2;
}

template<typename T>
  requires quantity<T>
auto Hilbert<T, T>::InnerProduct(T const& t1, T const& t2) -> InnerProductType {
  return t1 * t2;
}

template<typename T>
  requires quantity<T>
auto Hilbert<T, T>::Norm²(T const& t) -> Norm²Type {
  return t * t;
}

template<typename T>
  requires quantity<T>
auto Hilbert<T, T>::Norm(T const& t) -> NormType {
  return Abs(t);
}

template<typename T1, typename T2>
  requires hilbert<T1, T2>
auto Hilbert<T1, T2>::InnerProduct(T1 const& t1, T2 const& t2)
    -> InnerProductType {
  // Is there a better way to avoid recursion than to put our fingers inside
  // grassmann?
  return _grassmann::internal::InnerProduct(t1, t2);
}

template<typename T>
  requires hilbert<T, T>
auto Hilbert<T, T>::InnerProduct(T const& t1, T const& t2) -> InnerProductType {
  // Is there a better way to avoid recursion than to put our fingers inside
  // grassmann?
  return _grassmann::internal::InnerProduct(t1, t2);
}

template<typename T>
  requires hilbert<T, T>
auto Hilbert<T, T>::Norm²(T const& t) -> Norm²Type {
  return t.Norm²();
}

template<typename T>
  requires hilbert<T, T>
auto Hilbert<T, T>::Norm(T const& t) -> NormType {
  return t.Norm();
}

}  // namespace internal
}  // namespace _hilbert
}  // namespace geometry
}  // namespace principia
