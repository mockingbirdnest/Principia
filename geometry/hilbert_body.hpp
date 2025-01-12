#pragma once

#include "geometry/hilbert.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _hilbert {
namespace internal {

using namespace principia::quantities::_elementary_functions;

template<typename T1, typename T2>
  requires convertible_to_quantity<T1> && convertible_to_quantity<T2>
auto Hilbert<T1, T2>::InnerProduct(T1 const& t1, T2 const& t2)
    -> InnerProductType {
  return t1 * t2;
}

template<typename T>
  requires convertible_to_quantity<T>
auto Hilbert<T, T>::InnerProduct(T const& t1, T const& t2) -> InnerProductType {
  return t1 * t2;
}

template<typename T>
  requires convertible_to_quantity<T>
auto Hilbert<T, T>::Norm²(T const& t) -> Norm²Type {
  return t * t;
}

template<typename T>
  requires convertible_to_quantity<T>
auto Hilbert<T, T>::Norm(T const& t) -> NormType {
  return Abs(t);
}

#if !(_MSC_FULL_VER == 193'431'937 || \
      _MSC_FULL_VER == 193'431'942 || \
      _MSC_FULL_VER == 193'431'944 || \
      _MSC_FULL_VER == 193'532'216 || \
      _MSC_FULL_VER == 193'532'217 || \
      _MSC_FULL_VER == 193'632'532 || \
      _MSC_FULL_VER == 193'632'535 || \
      _MSC_FULL_VER == 193'732'822 || \
      _MSC_FULL_VER == 193'833'135 || \
      _MSC_FULL_VER == 193'933'523 || \
      _MSC_FULL_VER == 194'033'813 || \
      _MSC_FULL_VER == 194'134'120 || \
      _MSC_FULL_VER == 194'134'123 || \
      _MSC_FULL_VER == 194'234'435)
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
#endif

}  // namespace internal
}  // namespace _hilbert
}  // namespace geometry
}  // namespace principia
