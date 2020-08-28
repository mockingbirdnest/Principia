#pragma once

#include <type_traits>

#include "base/not_constructible.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/traits.hpp"

namespace principia {
namespace geometry {
namespace internal_hilbert {

using base::not_constructible;
using quantities::is_quantity;
using quantities::is_quantity_v;
using quantities::Product;
using quantities::Square;

// A trait that represents a Hilbert space, i.e., a space with an inner product
// and (possibly) a norm.  The struct Hilbert exports a type InnerProductType
// (the result of the inner product) and a function InnerProduct.  In addition,
// if only one parameter is given, or if the two parameters are identical, it
// also exports a type NormType (the result of the norm) and a function Norm.
template<typename T1, typename T2 = T1, typename = void>
struct Hilbert;

template<typename T1, typename T2>
struct Hilbert<T1, T2,
               std::void_t<std::enable_if_t<
                   std::conjunction_v<is_quantity<T1>, is_quantity<T2>>>>>
    : base::not_constructible {
  using InnerProductType = Product<T1, T2>;
  static InnerProductType InnerProduct(T1 const& t1, T2 const& t2);
};

template<typename T>
struct Hilbert<T, T, std::void_t<std::enable_if_t<is_quantity_v<T>>>>
    : base::not_constructible {
  using InnerProductType = Square<T>;
  static InnerProductType InnerProduct(T const& t1, T const& t2);

  using NormType = T;
  static NormType Norm(T const& t);
};

template<typename T1, typename T2>
struct Hilbert<T1, T2,
               std::void_t<decltype(InnerProduct(std::declval<T1>(),
                                                 std::declval<T2>()))>>
    : base::not_constructible {
  using InnerProductType =
      decltype(InnerProduct(std::declval<T1>(), std::declval<T2>()));
  static InnerProductType InnerProduct(T1 const& t1, T2 const& t2);
};

template<typename T>
struct Hilbert<T, T,
               std::void_t<decltype(InnerProduct(std::declval<T>(),
                                                 std::declval<T>()))>>
    : base::not_constructible {
  using InnerProductType =
      decltype(InnerProduct(std::declval<T>(), std::declval<T>()));
  static InnerProductType InnerProduct(T const& t1, T const& t2);

  using NormType = decltype(std::declval<T>().Norm());
  static NormType Norm(T const& t);
};

}  // namespace internal_hilbert

using internal_hilbert::Hilbert;

}  // namespace geometry
}  // namespace principia

#include "geometry/hilbert_body.hpp"
