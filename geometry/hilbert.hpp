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
using quantities::Product;
using quantities::SquareRoot;

template<typename T1, typename T2,
         typename = std::void_t<is_quantity<T1>, is_quantity<T2>>>
struct Hilbert : base::not_constructible {
  using InnerProductType = Product<T1, T2>;
  static InnerProductType InnerProduct(T1 const& t1, T2 const& t2);

  using NormType = std::conditional_t<std::is_same_v<T1, T2>, T1, void>;
  template<typename T = T1>
  static std::enable_if_t<std::is_same_v<T1, T2>, T> Norm(T const& t);
};

template<typename T1, typename T2>
struct Hilbert<T1, T2,
               std::void_t<decltype(InnerProduct(std::declval<T1>(),
                                                 std::declval<T2>()))>>
    : base::not_constructible {
  using InnerProductType =
      decltype(InnerProduct(std::declval<T1>(), std::declval<T2>()));
  static InnerProductType InnerProduct(T1 const& t1, T2 const& t2);

  using NormType = std::conditional_t<std::is_same_v<T1, T2>,
                                      decltype(std::declval<T1>().Norm()),
                                      void>;
  template<typename T = T1>
  static std::enable_if_t<std::is_same_v<T1, T2>,
                          decltype(std::declval<T>().Norm())> Norm(T const& t);
};

}  // namespace internal_hilbert

using internal_hilbert::Hilbert;

}  // namespace geometry
}  // namespace principia

#include "geometry/hilbert_body.hpp"
